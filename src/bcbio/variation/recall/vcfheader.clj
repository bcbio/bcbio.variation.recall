(ns bcbio.variation.recall.vcfheader
  "Create VCF headers using algorithms contained in Picard/Tribble tools.
   This does the best job of cleanly merging and organizing headers from
   multiple variant calling approaches."
  (:import [net.sf.picard.reference ReferenceSequenceFileFactory]
           [net.sf.picard.sam CreateSequenceDictionary]
           [org.broad.tribble AbstractFeatureReader]
           [org.broadinstitute.variant.vcf VCFCodec VCFUtils VCFHeader]
           [org.broadinstitute.variant.variantcontext.writer VariantContextWriterFactory
            Options])
  (:require [clojure.java.io :as io]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]))

(defn from-file
  "Retrieve header from input VCF file."
  [vcf-file]
  (println vcf-file)
  (with-open [vcf-reader (AbstractFeatureReader/getFeatureReader vcf-file (VCFCodec.) false)]
    (.getHeader vcf-reader)))

(defn create-ref-dict
  "Create reference dictionaries required by GATK and Picard.
   Requires samtools command to create *.fai if missing, since
   code to create these is no longer present in GATK."
  [ref-file]
  (let [dict-file (str (fsp/file-root ref-file) ".dict")]
    (when (itx/needs-run? dict-file)
      (.instanceMain (CreateSequenceDictionary.)
                     (into-array [(str "r=" ref-file) (str "o=" dict-file)])))
    dict-file))

(defn get-seq-dict
  "Retrieve Picard sequence dictionary from FASTA reference file."
  [ref-file]
  (create-ref-dict ref-file)
  (-> (io/file ref-file)
      ReferenceSequenceFileFactory/getReferenceSequenceFile
      .getSequenceDictionary))

(defn merge-from-files
  "Creates a merged VCF header from the supplied input VCFs."
  [orig-files ref-file out-file]
  (let [header-file (str (fsp/file-root out-file) "-header.vcf")
        headers (map from-file orig-files)]
    (with-open [vcf-writer (VariantContextWriterFactory/create (io/file header-file) (get-seq-dict ref-file)
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
