(ns bcbio.variation.recall.merge
  "Merge multiple VCF files together, running merge in parallel over genomic regions."
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.run.parallel :refer [rmap]]
            [bcbio.align.bed :as bed]
            [bcbio.variation.ensemble.prep :as eprep]
            [bcbio.variation.recall.split :as rsplit]
            [bcbio.variation.recall.vcfutils :as vcfutils]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [me.raynes.fs :as fs]))

(defn- region-merge
  "Perform merge of multiple sample VCFs within a given region."
  [vcf-files region work-dir final-file]
  (let [out-file (fsp/add-file-part final-file (eprep/region->safestr region)
                                    work-dir)
        vcf-file-str (string/join " " vcf-files)]
    (itx/run-cmd out-file
                 "bcftools merge -o ~{(vcfutils/bcftools-out-type out-file)} "
                 "-r ~{(eprep/region->samstr region)} ~{vcf-file-str} "
                 "> ~{out-file}")))

(defmulti concatenate-vcfs
  "Concatenate VCF files in supplied order, handling bgzip and plain text."
  (fn [_ out-file]
    (if (.endsWith out-file ".gz") :bgzip :default)))

(defmethod concatenate-vcfs :bgzip
  [vcf-files out-file]
  (eprep/bgzip-index-vcf (concatenate-vcfs vcf-files (subs out-file 0 (.lastIndexOf out-file ".gz")))
                         :remove-orig? true))

(defmethod concatenate-vcfs :default
  [vcf-files out-file]
  (let [str-vcf-files (string/join " " (rmap eprep/bgzip-index-vcf vcf-files))]
    (itx/run-cmd out-file
                 "vcfcat ~{str-vcf-files} > ~{out-file}")))

(defn combine-vcfs
  "Merge multiple VCF files together in parallel over genomic regions."
  [orig-vcf-files ref-file work-dir]
  (vcfutils/ensure-no-dup-samples orig-vcf-files)
  (let [vcf-files (rmap eprep/bgzip-index-vcf orig-vcf-files)
        region-bed (rsplit/group-pregions vcf-files ref-file work-dir)
        merge-dir (fsp/safe-mkdir (io/file (fs/parent region-bed) "merge"))
        out-file (fsp/add-file-part (fsp/remove-file-part region-bed "pregions")
                                    "merge" nil ".vcf.gz")
        merge-parts (reduce (fn [coll region]
                              (conj coll (region-merge vcf-files region merge-dir out-file)))
                            [] (bed/reader region-bed))]
    (concatenate-vcfs merge-parts out-file)))
