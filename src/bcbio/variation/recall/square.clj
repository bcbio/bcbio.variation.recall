(ns bcbio.variation.recall.square
  "Performing squaring off of variant call sets, recalling at all sample positions.
   This converts a merged dataset with no calls at positions not assessed in the
   sample, into a fully 'square' merged callset with reference calls at positions
   without evidence for a variant, distinguishing true no-calls from reference
   calls."
  (:require [bcbio.align.bam :as bam]
            [bcbio.align.cram :as cram]
            [bcbio.align.greads :as greads]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.variation.ensemble.prep :as eprep]
            [bcbio.variation.recall.clhelp :as clhelp]
            [bcbio.variation.recall.merge :as merge]
            [bcbio.variation.recall.vcfheader :as vcfheader]
            [bcbio.variation.recall.vcfutils :as vcfutils]
            [clojure.java.io :as io]
            [clojure.java.shell :refer [sh]]
            [clojure.string :as string]
            [clojure.tools.cli :refer [parse-opts]]
            [me.raynes.fs :as fs]
            [version-clj.core :refer [version-compare]]))

(defn subset-sample-region
  "Subset the input file to the given region and sample."
  [vcf-file sample region out-file]
  (itx/run-cmd out-file
               "bcftools view -O ~{(vcfutils/bcftools-out-type out-file)} "
               "-r ~{(eprep/region->samstr region)} -s ~{sample} "
               "~{(eprep/bgzip-index-vcf vcf-file)} > ~{out-file}")
  (eprep/bgzip-index-vcf out-file :remove-orig? true))

(defn- intersect-variants
  "Retrieve VCF variants present in both in-file and cmp-file."
  [in-file cmp-file ref-file out-file]
  (itx/run-cmd out-file
               "vcfintersect -r ~{ref-file} -i ~{cmp-file} ~{in-file} | "
               "bgzip > ~{out-file}")
  (eprep/bgzip-index-vcf out-file :remove-orig? true))

(defn- unique-variants
  "Retrieve variants from in-file not present in cmp-file."
  [in-file cmp-file ref-file out-file]
  (itx/run-cmd out-file
               "vcfintersect -v -r ~{ref-file} -i ~{cmp-file} ~{in-file} | "
               "bgzip > ~{out-file}")
  (eprep/bgzip-index-vcf out-file :remove-orig? true))

(defmulti recall-variants
  "Recall variants only at positions in provided input VCF file, using multiple callers."
  (fn [& args]
    (keyword (get (last args) :caller :freebayes))))

(defmethod recall-variants :freebayes
  ^{:doc "Perform variant recalling at specified positions with FreeBayes."}
  [sample region vcf-file bam-file ref-file out-file config]
  (let [sample-file (str (fsp/file-root out-file) "-samples.txt")]
    (spit sample-file sample)
    (itx/run-cmd out-file
                 "freebayes -b ~{bam-file} --variant-input ~{vcf-file} --only-use-input-alleles "
                 "--min-repeat-entropy 1 --experimental-gls "
                 " -f ~{ref-file} -r ~{(eprep/region->freebayes region)} -s ~{sample-file}  | "
                 "bgzip > ~{out-file}")
    (eprep/bgzip-index-vcf out-file :remove-orig? true)))

(defmulti platypus-filter
  "Perform post-filtration of platypus variant calls.
   Removes hard Q20 filter and replaces with NA12878/GiaB tuned depth
   and quality based filter. bgzips final output."
  (fn [approach f]
    (keyword approach)))

(defmethod platypus-filter :bcftools
  [_ orig-file]
  (let [out-file (str orig-file ".gz")]
    (itx/run-cmd out-file
                 " bcftools filter ~{orig-file} -s LOWDPQUAL -m '+' "
                 "-e '(((FR<=0.5)&&(TC<4)&&(%QUAL<20))||((TC<13)&&(%QUAL<10)))||((FR>0.5)&&((TC<4)&&(%QUAL<20)))' "
                 "| sed 's/\\tQ20\\t/\\tPASS\\t/' | bgzip -c > ~{out-file}")))

(defmethod platypus-filter :vcflib
  [_ orig-file]
  (let [out-file (str orig-file ".gz")]
    (itx/run-cmd out-file
                 "vcffilter ~{orig-file} -t LOWDPQUAL -A "
                 "--or -f 'FR < 0.55 & TC < 4 & QUAL < 20' -f 'FR < 0.55 & TC < 13 & QUAL < 10' "
                 "-f 'FR > 0.54 & TC < 4 & QUAL < 20'"
                 "| sed 's/\\tQ20\\t/\\tPASS\\t/' | sed 's/PASS,LOWDPQUAL/LOWDPQUAL/' "
                 "| bgzip -c > ~{out-file}")))

(defmethod recall-variants :platypus
  ^{:doc "Perform variant recalling at specified positions with Platypus."}
  [sample region vcf-file bam-file ref-file out-file config]
  (let [raw-out-file (string/replace out-file ".gz" "")]
    (when (itx/needs-run? out-file)
      (itx/run-cmd raw-out-file
                   "platypus callVariants --bamFiles=~{bam-file} --regions=~{(eprep/region->samstr region)} "
                   "--hapScoreThreshold 10 --scThreshold 0.99 --filteredReadsFrac 0.9 "
                   "--refFile=~{ref-file} --source=~{vcf-file} --minPosterior=0 --getVariantsFromBAMs=0 "
                   "--logFileName /dev/null --verbosity=1 --output ~{raw-out-file}")
      (platypus-filter :vcflib raw-out-file))
    (eprep/bgzip-index-vcf raw-out-file :remove-orig? true)))

(defn union-variants
  "Create a union of VCF variants from two files."
  [f1 f2 region ref-file out-file]
  (when (itx/needs-run? out-file)
    (vcfheader/with-merged [header-file [f1 f2] ref-file out-file]
      (itx/run-cmd out-file
                   "cat <(grep ^## ~{header-file}) <(vcfcombine ~{f1} ~{f2} | grep -v ^##) | "
                   "bgzip -c > ~{out-file}")))
  (eprep/bgzip-index-vcf out-file :remove-orig? true))

(defn- sample-by-region
  "Square off a specific sample in a genomic region, given all possible variants.
    - Subset to the current variant region.
    - Identify missing uncalled variants: create files of existing and missing variants.
    - Recall at missing positions with FreeBayes.
    - Merge original and recalled variants."
  [sample vcf-file bam-file union-vcf region ref-file out-file config]
  (let [work-dir (fsp/safe-mkdir (str (fsp/file-root out-file) "-work"))
        fnames (into {} (map (fn [x] [(keyword x) (str (io/file work-dir (format "%s.vcf.gz" x)))])
                             ["region" "existing" "needcall" "recall"]))]
    (when (itx/needs-run? out-file)
      (itx/with-temp-dir [tmp-dir (fs/parent out-file)]
        (let [region-bam-file (greads/subset-in-region bam-file ref-file region tmp-dir)]
          (subset-sample-region vcf-file sample region (:region fnames))
          (intersect-variants (:region fnames) union-vcf ref-file (:existing fnames))
          (unique-variants union-vcf (:region fnames) ref-file (:needcall fnames))
          (recall-variants sample region (:needcall fnames) region-bam-file ref-file (:recall fnames) config)
          (union-variants (:recall fnames) (:existing fnames) region ref-file out-file))))
    out-file))

(defn- sample-by-region-prep
  "Prepare for squaring off a sample in a region, setup out file and check conditions.
   We only can perform squaring off with a BAM file for the sample."
  [sample vcf-file bam-file union-vcf region ref-file out-dir config]
  (let [out-file (str (io/file (fsp/safe-mkdir (io/file out-dir (get region :chrom "nochrom")))
                               (format "%s-%s.vcf.gz" sample (eprep/region->safestr region))))]
    (cond
     (nil? bam-file) (subset-sample-region vcf-file sample region out-file)
     (itx/needs-run? out-file) (sample-by-region sample vcf-file bam-file union-vcf region ref-file out-file config)
     :else out-file)))

(defn by-region
  "Square off a genomic region, identifying variants from all samples and recalling at uncalled positions.
    - Identifies all called variants from all samples
    - For each sample, square off using `sample-by-region`
    - Merge all variant files in the region together."
  [vcf-files bam-files region ref-file dirs out-file config]
  (let [union-vcf (eprep/create-union vcf-files ref-file region (:union dirs))
        recall-vcfs (map (fn [[sample vcf-file]]
                           (sample-by-region-prep sample vcf-file (get bam-files sample)
                                                  union-vcf region ref-file (:square dirs) config))
                         (mapcat (fn [vf]
                                   (for [s (vcfutils/get-samples vf)]
                                     [s vf]))
                                 vcf-files))]
    (merge/region-merge :bcftools recall-vcfs region (:merge dirs) out-file)))

(defn sample-to-bam-map
  "Prepare a map of sample names to BAM files."
  [bam-files ref-file]
  (into {} (mapcat (fn [b]
                     (for [s (case (fs/extension b)
                               ".bam" (bam/sample-names b)
                               ".cram" (cram/sample-names b ref-file))]
                       [s b]))
                   bam-files)))

(defn- check-versions
  "Ensure we have up to date versions of required software for recalling."
  [config]
  (when (= :freebayes (:caller config))
    (let [version (-> (sh "freebayes")
                      :out
                      (string/split #"\n")
                      last
                      (string/split #":")
                      last
                      string/trim)
          want-version "v0.9.14-1"]
      (when (neg? (version-compare version want-version))
        (throw (Exception. (format "Require at least freebayes %s for recalling. Found %s"
                                   want-version version)))))))

(defn combine-vcfs
  "Combine VCF files with squaring off by recalling at uncalled variant positions."
  [orig-vcf-files bam-files ref-file out-file config]
  (let [dirs {:union (fsp/safe-mkdir (io/file (fs/parent out-file) "union"))
              :square (fsp/safe-mkdir (io/file (fs/parent out-file) "square"))}]
    (check-versions config)
    (merge/prep-by-region (fn [vcf-files region merge-dir]
                            (vcfutils/ensure-no-dup-samples vcf-files)
                            (by-region vcf-files (sample-to-bam-map bam-files ref-file)
                                       region ref-file (assoc dirs :merge merge-dir) out-file config))
                          orig-vcf-files ref-file out-file config)))

(defn- usage [options-summary]
  (->> ["Perform squaring off for a set of called VCF files, recalling at no-call positions in each sample."
        ""
        "Usage: bcbio-variation-recall square [options] out-file ref-file [<vcf, bam, cram, or list files>]"
        ""
        "  out-file:    VCF (or bgzipped VCF) file to write merged output to"
        "  ref-file:    FASTA format genome reference file"
        "  <remaining>: VCF files to recall and BAM or CRAM files for each sample. Can be specified "
        "               on the command line or as text files containing paths to files "
        "               for processing. VCFs can be single or multi-sample and BAM/CRAMs can be in "
        "               any order but each VCF sample must have an associated BAM/CRAM file to recall."
        ""
        "Options:"
        options-summary]
       (string/join \newline)))

(defn -main [& args]
  (let [caller-opts #{:freebayes :platypus}
        {:keys [options arguments errors summary]}
        (parse-opts args [["-c" "--cores CORES" "Number of cores to use" :default 1
                           :parse-fn #(Integer/parseInt %)]
                          ["-m" "--caller CALLER" (str "Calling method to use: "
                                                       (string/join ", " (map name caller-opts)))
                           :default "freebayes"
                           :parse-fn #(keyword %)
                           :validate [#(contains? caller-opts %)
                                      (str "Supported calling options: "
                                           (string/join ", " (map name caller-opts)))]]
                          ["-r" "--region REGION" "Genomic region to subset, in samtools format (chr1:100-200)"]
                          ["-h" "--help"]])]
    (cond
     (:help options) (clhelp/exit 0 (usage summary))
     errors (clhelp/exit 1 (clhelp/error-msg errors))
     (= 0 (count arguments)) (clhelp/exit 0 (usage summary))
     (< (count arguments) 3) (clhelp/exit 1 (usage summary))
     :else (let [[out-file ref-file & vcf-inputs] arguments
                 arg-files (clhelp/vcf-bam-args vcf-inputs)]
             (cond
              (not (empty? (:missing arg-files)))
              (clhelp/exit 1 (clhelp/error-msg (cons "Input files not found:" (:missing arg-files))))
              (or (not (fs/exists? ref-file)) (not (fs/file? ref-file)))
              (clhelp/exit 1 (clhelp/error-msg [(str "Reference file not found: " ref-file)]))
              :else
              (combine-vcfs (:vcf arg-files) (:bam arg-files) ref-file out-file options))))))
