(ns bcbio.variation.recall.split
  "Split variation calls into small genomic regions for parallel processing."
  (:require [bcbio.align.gref :as gref]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.variation.recall.vcfutils :as vcfutils]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [schema.core :as s]))

(s/defn ^:always-validate vcf-breakpoints
  "Prepare BED file of non-variant regions in the input VCF as parallel breakpoints.
   Uses bedtools to find covered regions by the VCF and subtracts this from the
   full reference genome to convert to non-covered/non-variant-call regions."
  ([vcf-file sample :- s/String ref-file work-dir]
     (let [merge-size 500
           split-dir (fsp/safe-mkdir (io/file work-dir "split" sample))
           sample-vcf (vcfutils/subset-to-sample vcf-file sample split-dir)
           out-file (fsp/add-file-part sample-vcf "splitpoints" split-dir ".bed")
           fai-file (gref/fasta-idx ref-file)]
       (itx/run-cmd out-file
                    "bedtools subtract "
                    "-a <(cut -f 1-2 ~{fai-file} | awk -F $'\t' '{OFS=FS} {print $1,0,$2}') "
                    "-b <(bedtools genomecov -i ~{sample-vcf} -g ~{fai-file} -bg | "
                    "     bedtools merge -d ~{merge-size}) "
                    "> ~{out-file}")))
  ([vcf-file ref-file work-dir]
     (vcf-breakpoints vcf-file (vcfutils/get-sample vcf-file) ref-file work-dir)))

(defn group-breakpoints
  "Prepare BED file of shared breakpoints for all samples in supplied files."
  [vcf-files ref-file work-dir]
  (let [split-dir (fsp/safe-mkdir (io/file work-dir "split"))
        vcf-samples (reduce (fn [coll vcf-file]
                              (concat coll (map (fn [s] [vcf-file s]) (vcfutils/get-samples vcf-file))))
                            [] vcf-files)
        out-file (fsp/add-file-part (first vcf-files) (format "combo-%s-splitpoints" (count vcf-samples))
                                    split-dir ".bed")
        bp-beds (map (fn [[v s]] (vcf-breakpoints v s ref-file work-dir)) vcf-samples)
        str-bp-beds (string/join " " bp-beds)]
    (itx/run-cmd out-file
                 "bedtools multiinter -i ~{str-bp-beds} | "
                 "awk -F $'\t' '{OFS=FS} $4 >= ~{(count bp-beds)} {print $1,$2,$3}' "
                 "> ~{out-file}")))
