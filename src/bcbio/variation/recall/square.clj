(ns bcbio.variation.recall.square
  "Performing squaring off of variant call sets, recalling at all sample positions.
   This converts a merged dataset with no calls at positions not assessed in the
   sample, into a fully 'square' merged callset with reference calls at positions
   without evidence for a variant, distinguishing true no-calls from reference
   calls."
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.variation.ensemble.prep :as eprep]
            [bcbio.variation.recall.merge :as merge]
            [clojure.java.io :as io]
            [me.raynes.fs :as fs]))

(defn by-region
  "Square off a genomic region, identifying variants from all samples and recalling at uncalled positions.
    - Identifies all called variants from all samples
    - For each sample, identify missing uncalled variants
    - Recall at these positions with FreeBayes
    - Merge original and recalled variants.
    - Merge all variant files in the region together."
  [vcf-files bam-files region ref-file out-dir union-dir out-file]
  (let [union-vcf (eprep/create-union vcf-files ref-file region union-dir)]
    (println union-vcf)
    union-vcf))

(defn combine-vcfs
  "Combine VCF files with squaring off by recalling at uncalled variant positions."
  [orig-vcf-files bam-files ref-file out-file config]
  (let [union-dir (fsp/safe-mkdir (io/file (fs/parent out-file) "union"))]
    (merge/prep-by-region (fn [vcf-files region out-dir]
                            (by-region vcf-files bam-files region ref-file out-dir union-dir out-file))
                          orig-vcf-files ref-file out-file config)))
