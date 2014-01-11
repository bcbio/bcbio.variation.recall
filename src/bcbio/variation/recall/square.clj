(ns bcbio.variation.recall.square
  "Performing squaring off of variant call sets, recalling at all sample positions.
   This converts a merged dataset with no calls at positions not assessed in the
   sample, into a fully 'square' merged callset with reference calls at positions
   without evidence for a variant, distinguishing true no-calls from reference
   calls."
  (:require [bcbio.align.bam :as bam]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.variation.ensemble.prep :as eprep]
            [bcbio.variation.recall.clhelp :as clhelp]
            [bcbio.variation.recall.merge :as merge]
            [bcbio.variation.recall.vcfutils :as vcfutils]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [clojure.tools.cli :refer [parse-opts]]
            [me.raynes.fs :as fs]))

(defn subset-sample-region
  "Subset the input file to the given region and sample."
  [vcf-file sample region out-file]
  (itx/run-cmd out-file
               "bcftools subset -o ~{(vcfutils/bcftools-out-type out-file)} "
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

(defn recall-variants
  "Recall variants only at positions in the provided input VCF file."
  [sample region vcf-file bam-file ref-file out-file]
  (let [raw-out-file (string/replace out-file ".gz" "")
        sample-file (str (fsp/file-root out-file) "-samples.txt")]
    (spit sample-file sample)
    (itx/run-cmd out-file
                 "freebayes -b ~{bam-file} -@ ~{vcf-file} -l -f ~{ref-file} "
                 "-r ~{(eprep/region->freebayes region)} -s ~{sample-file}  | "
                 "bgzip > ~{out-file}")
    (eprep/bgzip-index-vcf out-file :remove-orig? true)))

(defn union-variants
  "Create a union of VCF variants from two files."
  [f1 f2 region ref-file out-file]
  (itx/run-cmd out-file
               "vcfcombine ~{f1} ~{f2} | vcfcreatemulti | bgzip > ~{out-file}")
  (eprep/bgzip-index-vcf out-file :remove-orig? true))

(defn- sample-by-region
  "Square off a specific sample in a genomic region, given all possible variants.
    - Subset to the current variant region.
    - Identify missing uncalled variants: create files of existing and missing variants.
    - Recall at missing positions with FreeBayes.
    - Merge original and recalled variants."
  [sample vcf-file bam-file union-vcf region ref-file out-file]
  (let [work-dir (fsp/safe-mkdir (str (fsp/file-root out-file) "-work"))
        fnames (into {} (map (fn [x] [(keyword x) (str (io/file work-dir (format "%s.vcf.gz" x)))])
                                 ["region" "existing" "needcall" "recall"]))]
    (when (itx/needs-run? out-file)
      (subset-sample-region vcf-file sample region (:region fnames))
      (intersect-variants (:region fnames) union-vcf ref-file (:existing fnames))
      (unique-variants union-vcf (:region fnames) ref-file (:needcall fnames))
      (if (vcfutils/has-variants? (:needcall fnames))
        (do
          (recall-variants sample region (:needcall fnames) bam-file ref-file (:recall fnames))
          (union-variants (:existing fnames) (:recall fnames) region ref-file out-file))
        (merge/move-vcf (:region fnames) out-file)))
    out-file))

(defn- sample-by-region-prep
  "Prepare for squaring off a sample in a region, setup out file and check conditions.
   We only can perform squaring off with a BAM file for the sample."
  [sample vcf-file bam-file union-vcf region ref-file out-dir]
  (let [out-file (str (io/file (fsp/safe-mkdir (io/file out-dir (get region :chrom "nochrom")))
                               (format "%s-%s.vcf.gz" sample (eprep/region->safestr region))))]
    (cond
     (nil? bam-file) (subset-sample-region vcf-file sample region out-file)
     (itx/needs-run? out-file) (sample-by-region sample vcf-file bam-file union-vcf region ref-file out-file)
     :else out-file)))

(defn by-region
  "Square off a genomic region, identifying variants from all samples and recalling at uncalled positions.
    - Identifies all called variants from all samples
    - For each sample, square off using `sample-by-region`
    - Merge all variant files in the region together."
  [vcf-files bam-files region ref-file dirs out-file]
  (let [union-vcf (eprep/create-union vcf-files ref-file region (:union dirs))
        recall-vcfs (map (fn [[sample vcf-file]]
                           (sample-by-region-prep sample vcf-file (get bam-files sample)
                                                  union-vcf region ref-file (:square dirs)))
                         (mapcat (fn [vf]
                                   (for [s (vcfutils/get-samples vf)]
                                     [s vf]))
                                 vcf-files))]
    (merge/region-merge :bcftools recall-vcfs region (:merge dirs) out-file)))

(defn- sample-to-bam-map
  "Prepare a map of sample names to BAM files."
  [bam-files]
  (into {} (mapcat (fn [b]
                     (for [s (bam/sample-names b)]
                       [s b]))
                   bam-files)))

(defn combine-vcfs
  "Combine VCF files with squaring off by recalling at uncalled variant positions."
  [orig-vcf-files bam-files ref-file out-file config]
  (let [dirs {:union (fsp/safe-mkdir (io/file (fs/parent out-file) "union"))
              :square (fsp/safe-mkdir (io/file (fs/parent out-file) "square"))}]
    (merge/prep-by-region (fn [vcf-files region merge-dir]
                            (by-region vcf-files (sample-to-bam-map bam-files)
                                       region ref-file (assoc dirs :merge merge-dir) out-file))
                          orig-vcf-files ref-file out-file config)))

(defn- usage [options-summary]
  (->> ["Perform squaring off for a set of called VCF files, recalling at no-call positions in each sample."
        ""
        "Usage: bcbio-variation-recall square [options] out-file ref-file [<vcf-files, bam-files, or list-files>]"
        ""
        "  out-file:    VCF (or bgzipped VCF) file to write merged output to"
        "  ref-file:    FASTA format genome reference file"
        "  <remaining>: VCF files to recall and BAM files for each sample. Can be specified "
        "               on the command line or as text files containing paths to files "
        "               for processing. VCFs can be single or multi-sample and BAMs can be in "
        "               any order but each VCF sample must have an associated BAM file to recall."
        ""
        "Options:"
        options-summary]
       (string/join \newline)))

(defn -main [& args]
  (let [{:keys [options arguments errors summary]}
        (parse-opts args [["-c" "--cores CORES" "Number of cores to use" :default 1
                           :parse-fn #(Integer/parseInt %)]
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
