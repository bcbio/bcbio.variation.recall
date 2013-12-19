(ns bcbio.variation.recall.vcfutils
  "Utilities for manipulating VCF files"
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [schema.core :as s]))

(defn pog-reader
  "Plain or gzip input reader."
  [f]
  (if (.endsWith f ".gz")
    (io/reader (java.util.zip.GZIPInputStream. (io/input-stream f)))
    (io/reader f)))

(defn get-samples
  "Retrieve sample names from a VCF file."
  [vcf-file]
  (with-open [rdr (pog-reader vcf-file)]
    (let [line (first (drop-while #(not (.startsWith % "#CHROM")) (line-seq rdr)))]
      (vec (drop 9 (string/split (string/trimr line) #"\t"))))))

(defn get-sample
  "Retrieve VCF sample name from a single sample input file."
  [vcf-file]
  (let [samples (get-samples vcf-file)]
    (when (= 1 (count samples))
      (first samples))))

(defn duplicate-samples
  "Retrieve any duplicate samples in the input VCF files."
  [vcf-files]
  (->> vcf-files
       (mapcat get-samples)
       frequencies
       (filter (fn [[x n]] (> n 1)))
       (map first)))

(defn ensure-no-dup-samples
  "Ensure there are sample name duplicates in the VCF files"
  [vcf-files]
  (let [dups (duplicate-samples vcf-files)]
    (assert (empty? dups) (format "Found duplicate samples %s in VCFs for merging: %s"
                                  (vec dups) (vec vcf-files)))))

(defn bcftools-out-type
  [f]
  (if (.endsWith f ".gz") "z" "v"))

(defn subset-to-sample*
  "Do the actual work of subsetting a file, assumes multisample VCF"
  [vcf-file sample out-file]
  (let [out-type (bcftools-out-type out-file)]
    (itx/run-cmd out-file
                 "bcftools subset -s ~{sample} -o ~{out-type} ~{vcf-file} > ~{out-file}"))
  out-file)

(defn subset-to-sample
  "Subset a VCF file, retrieving a single sample from multisample VCF.
   If we already have a single sample VCF, return that."
  [vcf-file sample out-dir]
  (if (> (count (get-samples vcf-file)) 1)
    (subset-to-sample* vcf-file sample (fsp/add-file-part vcf-file sample out-dir))
    vcf-file))
