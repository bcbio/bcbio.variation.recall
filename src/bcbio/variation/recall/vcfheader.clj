(ns bcbio.variation.recall.vcfheader
  "Create VCF headers using algorithms contained in Picard/Tribble tools.
   This does the best job of cleanly merging and organizing headers from
   multiple variant calling approaches."
  (:import [htsjdk.variant.vcf VCFUtils VCFHeader]
           [htsjdk.variant.variantcontext.writer VariantContextWriterFactory])
  (:require [clojure.java.io :as io]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.variation.variantcontext :as gvc]))

(defn merge-from-files
  "Creates a merged VCF header from the supplied input VCFs."
  [orig-files ref-file out-file]
  (let [header-file (str (fsp/file-root out-file) "-header.vcf")
        headers (map gvc/get-vcf-header orig-files)]
    (with-open [vcf-writer (VariantContextWriterFactory/create (io/file header-file) (gvc/get-seq-dict ref-file)
                                                               VariantContextWriterFactory/NO_OPTIONS)]
      (.writeHeader vcf-writer (VCFHeader. (VCFUtils/smartMergeHeaders headers false))))
    header-file))

(defmacro with-merged
  "Create a merged VCF header file from input VCFs, deleted on completion."
  [[header-file orig-files ref-file out-file] & body]
  `(let [~header-file (merge-from-files ~orig-files ~ref-file ~out-file)]
     (try
       (let [out# (do ~@body)]
         out#)
       (finally
        (fsp/remove-path ~header-file)))))
