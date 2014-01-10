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
            [bcbio.variation.recall.merge :as merge]
            [bcbio.variation.recall.vcfutils :as vcfutils]
            [clojure.java.io :as io]
            [me.raynes.fs :as fs]))

(defn subset-sample-region
  "Subset the input file to the given region and sample."
  [vcf-file sample region out-file]
  (itx/run-cmd out-file
               "bcftools subset -o ~{(vcfutils/bcftools-out-type out-file)} "
               "-r ~{(eprep/region->samstr region)} -s ~{sample} "
               "~{(eprep/bgzip-index-vcf vcf-file)} > ~{out-file}")
  (eprep/bgzip-index-vcf out-file :remove-orig? true))

(defn- sample-by-region
  "Square off a specific sample in a genomic region, given all possible variants.
    - Subset to the current variant region.
    - Identify missing uncalled variants: create files of existing and missing variants.
    - Recall at missing positions with FreeBayes.
    - Merge original and recalled variants."
  [sample vcf-file bam-file union-vcf region ref-file out-file]
  (println out-file))

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
    (vec recall-vcfs)
    ;(merge/region-merge :bcftools recall-vcfs region (:merge dirs) out-file)
    union-vcf))

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
