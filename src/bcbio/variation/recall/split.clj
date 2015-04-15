(ns bcbio.variation.recall.split
  "Split variation calls into small genomic regions for parallel processing."
  (:require [bcbio.align.gref :as gref]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.run.parallel :refer [rmap]]
            [bcbio.variation.ensemble.prep :as eprep]
            [bcbio.variation.recall.vcfutils :as vcfutils]
            [clojure.core.strint :refer [<<]]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [ordered.map :refer [ordered-map]]
            [me.raynes.fs :as fs]))

(def ^{:private true} merge-size 25000)

(defn- region->bed
  "Convert specified input region into a BED file"
  [region out-file]
  (when region
    (if (fs/file? region)
      (itx/safe-copy (io/file region) (io/file out-file))
      (spit out-file (str (string/join "\t" (string/split region #"[:-]")) "\n")))
    out-file))

(defn- fai->bed
  "Convert fasta reference file to genome, potentially subsetting by the region of interest."
  [fai-file region base-file]
  (let [region-file (region->bed region (fsp/add-file-part base-file "inputregion" nil ".bed"))]
    (str (<< "cut -f 1-2 ~{fai-file} | awk -F $'\\t' '{OFS=FS} {print $1,0,$2}'")
         (if region-file (<< " | bedtools intersect -a stdin -b ~{region-file} | sort -k1,1 -k2,2n") ""))))

(defn- regional-faibed
  "Convert a full genome fai file into the current genome region in BED format."
  [fai-file region base-file]
  (let [out-file (fsp/add-file-part base-file "faibed" nil ".bed")]
    (itx/run-cmd out-file
                 "~{(fai->bed fai-file region base-file)} > ~{out-file}")))

(defn- vcf-breakpoints
  "Prepare BED file of non-variant regions in the input VCF as parallel breakpoints.
   Uses bedtools to find covered regions by the VCF and subtracts this from the
   full reference genome to convert to non-covered/non-variant-call regions."
  ([vcf-file sample ref-file split-dir work-dir config]
     (let [split-dir (fsp/safe-mkdir (io/file split-dir sample))
           sample-vcf (eprep/bgzip-index-vcf (vcfutils/subset-to-sample vcf-file sample split-dir (:region config)))
           out-file (fsp/add-file-part sample-vcf "splitpoints" split-dir ".bed")
           fai-file (gref/fasta-idx ref-file)
           fai-bed-file (regional-faibed fai-file (:region config) out-file)]
       {:vcf sample-vcf
        :bed (if (vcfutils/has-variants? sample-vcf)
               (itx/run-cmd out-file
                            "bedtools subtract "
                            "-a ~{fai-bed-file} "
                            "-b <(bedtools slop -header -i ~{sample-vcf} -b 10 -g ~{fai-file}"
                            "     | bedtools merge -d ~{merge-size} -i stdin)"
                            " | sort -k1,1 -k2,2n "
                            "> ~{out-file}")
               (region->bed (:region config) out-file))}))
  ([vcf-file ref-file split-dir work-dir config]
     (vcf-breakpoints vcf-file (vcfutils/get-sample vcf-file) ref-file split-dir work-dir config)))

(defn- group-breakpoints
  "Prepare BED file of shared breakpoints for all samples in supplied files."
  [vcf-files ref-file work-dir config]
  (let [split-dir (fsp/safe-mkdir (io/file work-dir "split" (vcfutils/region->fileext (:region config))))
        vcf-samples (reduce (fn [coll vcf-file]
                              (concat coll (map (fn [s] [vcf-file s]) (vcfutils/get-samples vcf-file))))
                            [] vcf-files)
        out-file (fsp/add-file-part (first vcf-files) (format "combo-%s-splitpoints" (count vcf-samples))
                                    split-dir ".bed")
        bp-info (rmap (fn [[v s]] (vcf-breakpoints v s ref-file work-dir split-dir config)) vcf-samples
                      (:cores config))
        str-bp-beds (when (every? #(> (fs/size %) 0) (map :bed bp-info))
                          (string/join " " (map :bed bp-info)))]
    [(map :vcf bp-info)
     (if str-bp-beds
       (itx/run-cmd out-file
                    "bedtools multiinter -i ~{str-bp-beds} | "
                    "awk -F $'\\t' '{OFS=FS} $4 >= ~{(count bp-info)} {print $1,$2,$3}' "
                    "> ~{out-file}")
       out-file)]))

(defn group-pregions
  "Provide BED file of common analysis regions across all VCFs for parallel analysis.
   Creates pad 200bp on either side of variants to pull in surrounding region for re-analysis."
  [orig-vcf-files ref-file work-dir config]
  (let [pad-size 200
        [vcf-files bp-file] (group-breakpoints orig-vcf-files ref-file work-dir config)
        out-file (fsp/add-file-part (fsp/remove-file-part bp-file "splitpoints") "pregions")
        fai-file (gref/fasta-idx ref-file)]
    [vcf-files
     (if (and (fs/exists? bp-file) (> (fs/size bp-file) 0))
       (itx/run-cmd out-file
                    "bedtools subtract -a <(~{(fai->bed fai-file (:region config) out-file)}) -b ~{bp-file} | "
                    "bedtools slop -b ~{pad-size} -g ~{fai-file} -i | "
                    "bedtools merge -d ~{merge-size} -i stdin "
                    "> ~{out-file}")
       (regional-faibed fai-file (:region config) out-file))]))
