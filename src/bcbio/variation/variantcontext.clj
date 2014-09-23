(ns bcbio.variation.variantcontext
  "Helper functions to retrieve information from htsjdk VariantContext
   objects, which represent variant data stored in VCF files."
  (:import [htsjdk.samtools.reference ReferenceSequenceFileFactory]
           [htsjdk.tribble AbstractFeatureReader]
           [htsjdk.variant.vcf
            VCFCodec VCFUtils VCFHeader VCFFilterHeaderLine]
           [htsjdk.variant.variantcontext VariantContextBuilder
            GenotypeBuilder GenotypesContext]
           [htsjdk.variant.variantcontext.writer VariantContextWriterBuilder
            Options]
           [picard.sam CreateSequenceDictionary]
           [java.util EnumSet])
  (:use [clojure.set :only [intersection union]]
        [lazymap.core :only [lazy-hash-map]]
        [ordered.set :only [ordered-set]])
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.variation.ensemble.prep :as eprep]))

;; ## Represent VariantContext objects
;;
;; Provide simple map-based access to important attributes of
;; VariantContexts. There are 3 useful levels of abstraction:
;;
;;  - VariantContext: Details about a variation. This captures a
;;    single line in a VCF file
;;  - Genotype: An individual genotype for a sample, at a variant position.
;;  - Allele: The actual alleles at a genotype.

(defn from-genotype
  "Represent a sample genotype including alleles.
   :genotype stores the original java genotype object for direct access."
  [g]
  (lazy-hash-map
   :sample-name (.getSampleName g)
   :qual (.getPhredScaledQual g)
   :type (-> g .getType .name)
   :phased? (.isPhased g)
   :attributes (merge {"DP" (.getDP g) "AD" (vec (.getAD g))
                       "GQ" (.getGQ g) "PL" (vec (.getPL g))}
                      (into {} (.getExtendedAttributes g)))
   :alleles (vec (.getAlleles g))
   :genotype g))

(defn from-vc
  "Provide a top level map of information from a variant context.
   :vc stores the original java VariantContext object for direct access."
  [vc]
  (lazy-hash-map
   :chr (.getChr vc)
   :start (.getStart vc)
   :end (.getEnd vc)
   :id (when (.hasID vc) (.getID vc))
   :ref-allele (.getReference vc)
   :alt-alleles (vec (.getAlternateAlleles vc))
   :type (-> vc .getType .name)
   :filters (set (.getFilters vc))
   :attributes (into {} (.getAttributes vc))
   :qual (.getPhredScaledQual vc)
   :num-samples (.getNSamples vc)
   :genotypes (map from-genotype
                   (-> vc .getGenotypes .toArray vec))
   :vc vc))

;; ## Parsing VCF files

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

(defn get-vcf-source
  "Create a Tribble FeatureSource for VCF file, using bgzipped tabix indexed VCFs."
  ([in-file]
     (eprep/bgzip-index-vcf in-file)
     (AbstractFeatureReader/getFeatureReader in-file (VCFCodec.) false))
  ([in-file ref-file]
     (get-vcf-source in-file)))

(defn get-vcf-iterator
  "Create an iterator over VCF VariantContexts."
  ([in-file]
     (.iterator (get-vcf-source in-file)))
  ([in-file ref-file]
     (get-vcf-iterator in-file)))

(defn variants-in-region
  "Retrieve variants located in potentially multiple variant files"
  ([retriever vc]
     (variants-in-region retriever (:chr vc) (:start vc) (:end vc)))
  ([retriever space start end]
     (letfn [(get-vcs-in-source [[source fname]]
               (with-open [vcf-iter (.query source space start end)]
                 (doall (map #(assoc (from-vc %) :fname fname) (iterator-seq vcf-iter)))))]
       (mapcat get-vcs-in-source (map vector (:sources retriever) (:fnames retriever))))))

(defrecord VariantRetriever [sources fnames]
  java.io.Closeable
  (close [_]
    (doseq [x sources]
      (.close x))))

(defn get-vcf-retriever
  "Indexed variant file retrieval for zero to multiple files with clean handle closing."
  [ref & vcf-files]
  (let [fnames (remove nil? vcf-files)]
    (VariantRetriever. (map #(get-vcf-source % ref) fnames)
                       fnames)))

(defn parse-vcf
  "Lazy iterator of VariantContext information from VCF file."
  [vcf-source]
  (map from-vc (iterator-seq (.iterator vcf-source))))

;; ## Writing VCF files

(defn get-vcf-header
  "Retrieve header from input VCF file."
  [vcf-file]
  (with-open [vcf-reader (AbstractFeatureReader/getFeatureReader vcf-file (VCFCodec.) false)]
    (.getHeader vcf-reader)))

(defn merge-headers
  [& merge-files]
  (fn [_ header]
    (VCFHeader. (VCFUtils/smartMergeHeaders (cons header (map get-vcf-header merge-files))
                                            true)
                (.getGenotypeSamples header))))

(defn write-vcf-w-template
  "Write VCF output files starting with an original input template VCF.
   Handles writing to multiple VCF files simultaneously with the different
   file handles represented as keywords. This allows lazy splitting of VCF files:
   `vc-iter` is a lazy sequence of `(writer-keyword variant-context)`.
   `out-file-map` is a map of writer-keywords to output filenames."
  [tmpl-file out-file-map vc-iter & {:keys [header-update-fn]}]
  (letfn [(make-vcf-writer [f]
            (-> (VariantContextWriterBuilder. )
                (.setOutputFile f)
                (.setOption Options/INDEX_ON_THE_FLY)
                (.setOption Options/ALLOW_MISSING_FIELDS_IN_HEADER)
                .build))
          (convert-to-output [info]
            [(if (and (coll? info) (= 2 (count info))) (first info) :out)
             (if (coll? info) (last info) info)])]
    (itx/with-tx-files [tx-out-files out-file-map (keys out-file-map) [".idx" ".tbi"]]
      (let [tmpl-header (get-vcf-header tmpl-file)
            writer-map (zipmap (keys tx-out-files)
                               (map make-vcf-writer (vals tx-out-files)))]
        (doseq [[key out-vcf] writer-map]
          (.writeHeader out-vcf (if-not (nil? header-update-fn)
                                  (header-update-fn key tmpl-header)
                                  tmpl-header)))
        (doseq [[fkey item] (map convert-to-output vc-iter)]
          (let [ready-vc (if (and (map? item) (contains? item :vc)) (:vc item) item)]
            (when-not (nil? ready-vc)
              (.add (get writer-map fkey) ready-vc))))
        (doseq [x (vals writer-map)]
          (.close x))))))

;; ## Utilities

(defn remove-filter
  "Remove filter from a variant context."
  [vc]
  (-> (VariantContextBuilder. vc) .passFilters .make))
