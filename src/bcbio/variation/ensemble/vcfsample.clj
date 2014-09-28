(ns bcbio.variation.ensemble.vcfsample
  "Sort VCF sample columns to have a consistent order between multiple inputs.
   Variant callers order called outputs differently and this ensures they are
   consistent to feed into ensemble calling."
  (:import [htsjdk.samtools.util BlockCompressedInputStream BlockCompressedOutputStream])
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.variation.ensemble.prep :as eprep]
            [bcbio.variation.variantcontext :as vc]
            [clojure.java.io :as io]
            [clojure.string :as string]))

(defn- sort-sample-line
  "Sort samples in a VCF line using reordered indexes from calculate-reorder."
  [line reorder]
  (let [[keep samples] (split-at 9 (string/split line #"\t"))]
    (string/join "\t"
                 (concat keep
                         (->> samples
                              (map-indexed vector)
                              (sort-by (fn [[i x]] (get reorder i)))
                              (map second))))))

(defn- calculate-reorder
  "Create a dictionary of sample indexes in the original VCF to those desired in the output."
  [orig-order want-order]
  (let [want-indexes (reduce (fn [coll [i x]]
                               (assoc coll x i))
                             {} (map-indexed vector want-order))]
    (reduce (fn [coll [i x]]
              (assoc coll i (get want-indexes x)))
            {} (map-indexed vector orig-order))))

(defn unique-work-file
  "Create a work file with unique name in case of shared base names."
  [orig-file ext all-files work-dir]
  (let [cmp-files (remove #(= % orig-file) all-files)
        parts (reverse (string/split orig-file #"/"))
        unique-file (loop [i 1]
                      (let [cur-name (string/join "-" (reverse (take i parts)))]
                        (if (not-any? #(.endsWith % cur-name) cmp-files)
                          cur-name
                          (recur (inc i)))))]
    (fsp/add-file-part unique-file ext work-dir)))

(defn- sort-samples
  "Sort samples in a VCF file, moving from orig-order to want-order."
  [vcf-file orig-order want-order all-vcfs work-dir]
  (let [out-file (unique-work-file vcf-file "ssort" all-vcfs work-dir)
        sample-reorder (calculate-reorder orig-order want-order)]
    (with-open [rdr (io/reader (BlockCompressedInputStream. (io/file vcf-file)))
                wtr (io/writer (BlockCompressedOutputStream. (io/file out-file)))]
      (doseq [line (line-seq rdr)]
        (.write wtr (str (if (.startsWith line "##")
                           line
                           (sort-sample-line line sample-reorder))
                         "\n"))))
    (eprep/bgzip-index-vcf out-file)
    out-file))

(defn- maybe-sort-names
  "Extract sample names for the current file and do sorting if needed."
  [vcf-file sorder all-vcfs work-dir]
  (let [cur-sorder (vc/get-vcf-samples vcf-file)]
    (if (not= cur-sorder sorder)
      (sort-samples vcf-file cur-sorder sorder all-vcfs work-dir)
      vcf-file)))

(defn consistent-order
  "Ensure the set of VCF files have a consistent sample order relative to the first."
  [vcf-files work-dir]
  (fsp/safe-mkdir work-dir)
  (let [sorder (vc/get-vcf-samples (first vcf-files))]
    (cons (first vcf-files)
          (map #(maybe-sort-names % sorder vcf-files work-dir) (rest vcf-files)))))
