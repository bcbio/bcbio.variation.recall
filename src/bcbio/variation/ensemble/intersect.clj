(ns bcbio.variation.ensemble.intersect
  "Intersection based Ensemble calling approaches using inputs from multiple callers."
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.run.parallel :refer [rmap]]
            [bcbio.variation.ensemble.prep :as eprep]
            [bcbio.variation.recall.clhelp :as clhelp]
            [clojure.tools.cli :refer [parse-opts]]
            [clojure.string :as string]
            [me.raynes.fs :as fs]))

(defn- intersect-vcfs
  [vcf-files base-file options]
  (let [out-file (str (fsp/file-root base-file) "-intersect.txt")
        vcf-file-str (string/join " " vcf-files)
        numpass (get options :numpass 2)]
    (itx/run-cmd out-file
                 "bcftools isec -n '+~{numpass}' -o ~{out-file} ~{vcf-file-str}")))

(defn ensemble-vcfs
  "Calculate ensemble calls using intersection counting from a set of input VCFs"
  [orig-vcf-files ref-file out-file options]
  (fsp/safe-mkdir (fs/parent out-file))
  (let [vcf-files (rmap eprep/bgzip-index-vcf orig-vcf-files (:cores options))
        isec-file (intersect-vcfs vcf-files out-file options)]
    (println isec-file)
    out-file))

(defn- usage [options-summary]
  (->> ["Ensemble calling for samples: combine multiple VCF caller outputs into a single callset."
        ""
        "Usage: bcbio-variation-recall ensemble [options] out-file ref-file [<vcf-files, bam-files, or list-files>]"
        ""
        "   out-file:   VCF (or bgzipped VCF) file to write merged output to"
        "   ref-file:   FASTA format genome reference file"
        "  <remaining>: VCF files to include for building a final ensemble callset. "
        "               Specify on the command line or as text files containing paths to files "
        "               for processing. Variants are extracted for the final file in order supplied. "
        "               VCFs can be single or multi-sample. "
        ""
        "Options:"
        options-summary]
       (string/join \newline)))

(defn -main [& args]
  (let [{:keys [options arguments errors summary]}
        (parse-opts args [["-c" "--cores CORES" "Number of cores to use" :default 1
                           :parse-fn #(Integer/parseInt %)]
                          ["-n", "--numpass NUMPASS" "Number of callers a variant should be present in to pass"
                           :default 2 :parse-fn #(Integer/parseInt %)]
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
              (ensemble-vcfs (:vcf arg-files) ref-file out-file options))))))
