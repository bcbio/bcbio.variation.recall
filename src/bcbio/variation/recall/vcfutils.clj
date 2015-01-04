(ns bcbio.variation.recall.vcfutils
  "Utilities for manipulating VCF files"
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [me.raynes.fs :as fs]
            [taoensso.timbre :as timbre]))

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
      (if line
        (vec (drop 9 (string/split (string/trimr line) #"\t")))
        (throw (Exception. (format "Did not find #CHROM line in input VCF file %s" vcf-file)))))))

(defn get-sample
  "Retrieve VCF sample name from a single sample input file."
  [vcf-file]
  (let [samples (get-samples vcf-file)]
    (when (= 1 (count samples))
      (first samples))))

(defn has-variants?
  "Does a VCF file have any variants."
  [vcf-file]
  (with-open [rdr (pog-reader vcf-file)]
    (not (nil? (first (drop-while #(.startsWith % "#") (line-seq rdr)))))))

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
    (when-not (empty? dups)
      (let [e (Exception. (format "Found duplicate samples %s in VCFs for merging: %s"
                                  (vec dups) (vec vcf-files)))]
        (timbre/error e)
        (throw e)))))

(defn bcftools-out-type
  [f]
  (if (.endsWith f ".gz") "z" "v"))

(defn- region->bcftools
  "Convert a region specification to a bcftools option, handling chr1:1-100 and BED files"
  [region]
  (cond
   (and (fs/file? region) (fs/exists? region)) (str "-R " region)
   region (str "-r " region)
   :else ""))

(defn region->fileext
  "Convert a region to a stable file extension, handling chr1:1-100 and BED files"
  ([region-orig prefix]
     (let [region (if (and region-orig (fs/file? region-orig) (fs/exists? region-orig))
                    (with-open [rdr (io/reader region-orig)]
                      (->> (line-seq rdr)
                           (remove #(.startsWith % "track"))
                           (remove #(.startsWith % "#"))
                           first
                           (#(string/split % #"\t"))
                           (take 3)
                           (string/join "_")))
                    region-orig)]
       (if region
         (str prefix (string/replace region #"[-:]" "_"))
         "")))
  ([region-orig]
     (region->fileext region-orig "")))

(defn- subset-to-sample*
  "Do the actual work of subsetting a file, assumes multisample VCF"
  [vcf-file sample out-file region]
  (let [out-type (bcftools-out-type out-file)
        region_str (region->bcftools region)]
    (itx/run-cmd out-file
                 "bcftools view -s ~{sample} ~{region_str} -O ~{out-type} ~{vcf-file} > ~{out-file}"))
  out-file)

(defn subset-to-sample
  "Subset a VCF file, retrieving a single sample from multisample VCF.
   If we already have a single sample VCF, return that."
  ([vcf-file sample out-dir region]
     (if (or region (> (count (get-samples vcf-file)) 1))
       (let [ext (format "%s%s" sample (region->fileext region "-"))
             out-file (fsp/add-file-part vcf-file ext out-dir)]
         (subset-to-sample* vcf-file sample out-file region))
       vcf-file))
  ([vcf-file sample out-dir]
     (subset-to-sample vcf-file sample out-dir nil)))
