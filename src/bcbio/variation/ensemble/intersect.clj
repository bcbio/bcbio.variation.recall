(ns bcbio.variation.ensemble.intersect
  "Intersection based Ensemble calling approaches using inputs from multiple callers."
  (:import  [htsjdk.variant.vcf VCFInfoHeaderLine VCFHeaderLineType VCFHeaderLineCount])
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.run.clhelp :as clhelp]
            [bcbio.run.parallel :refer [rmap]]
            [bcbio.variation.ensemble.prep :as eprep]
            [bcbio.variation.ensemble.vcfsample :as vcfsample]
            [bcbio.variation.variantcontext :as vc]
            [clojure.java.io :as io]
            [clojure.tools.cli :refer [parse-opts]]
            [clojure.string :as string]
            [me.raynes.fs :as fs]))

(defn- intersect-vcfs
  [vcf-files work-dir base-file options]
  (let [out-file (str (fsp/add-file-part base-file "intersect" work-dir ".txt"))
        vcf-file-str (string/join " " vcf-files)
        numpass (get options :numpass 2)]
    (itx/run-cmd out-file
                 "bcftools isec -n '+~{numpass}' -o ~{out-file} ~{vcf-file-str}")))

(defn- parse-isec-line
  "Convert a line from bcftools isec into coordinates and the file to retrieve."
  [line]
  (let [[chrom str-start refa alta intersects] (string/split line #"\t")
        start (dec (Integer/parseInt str-start))]
    {:chr chrom :start start :refa (clojure.string/upper-case refa)
     :alta (map clojure.string/upper-case (string/split alta #","))
     :end (->> (cons refa (string/split alta #","))
              (map count)
              (apply max)
              (+ start))
     :vc-indices (->> intersects
                      (map-indexed (fn [i x] (when (= \1 x) i)))
                      (remove nil?))}))

(defn- ann-vc-w-names
  "Annotate a variant context with information contributing caller names."
  [context names]
  (if names
    (vc/vc-add-attr context "CALLERS" (string/join "," names))
    context))

(defn- get-rep-vc
  "Retrieve the best representative variant context for a passing ensemble variant."
  [vc-getter vcf-files names]
  (fn [line]
    (let [rep-file (nth vcf-files (first (:vc-indices line)))
          rep-vcs (->> (vc/variants-in-region vc-getter line)
                       (filter #(= rep-file (:fname %)))
                       (filter #(= (:start line) (dec (:start %))))
                       (filter #(= (:refa line) (:ref-allele %)))
                       (filter #(= (:alta line) (:alt-alleles %))))
          cur-names (when (seq names) (map (partial nth names) (:vc-indices line)))]
      (if (> (count rep-vcs) 0)
        (-> rep-vcs first :vc vc/remove-filter (ann-vc-w-names cur-names))
        (throw (Exception. (format "Problem retrieving reference variant for %s: %s" line
                                   (vec (map #(select-keys % [:chr :start :fname]) rep-vcs)))))))))

(defn- maybe-nofiltered
  "Potentially remove filtered variants from inputs into Ensemble calling.
   Enables inputs to Ensemble calling to use all variants or only
   those that pass an initial filtration step."
  [vcf-file all-files work-dir options]
  (let [out-file (eprep/unique-work-file vcf-file "nofilter" all-files work-dir)]
    (if (:nofiltered options)
      (do
        (itx/run-cmd out-file
                     "bcftools view -f 'PASS,.' ~{vcf-file} -O z -o ~{out-file}")
        (eprep/bgzip-index-vcf out-file))
      vcf-file)))

(defn ensemble-vcfs
  "Calculate ensemble calls using intersection counting from a set of input VCFs"
  [orig-vcf-files ref-file orig-out-file options]
  (let [out-file (if (.endsWith orig-out-file ".gz") orig-out-file (format "%s.gz" orig-out-file))]
    (fsp/safe-mkdir (fs/parent out-file))
    (when (itx/needs-run? out-file)
      (let [bg-vcf-files (rmap eprep/bgzip-index-vcf orig-vcf-files (:cores options))
            work-dir (fsp/safe-mkdir (str (fsp/file-root out-file) "-work"))
            vcf-files (vcfsample/consistent-order (rmap #(maybe-nofiltered % bg-vcf-files work-dir options)
                                                        bg-vcf-files)
                                                  work-dir)
            isec-file (intersect-vcfs vcf-files work-dir out-file options)]
        (with-open [rdr (io/reader isec-file)
                    vc-getter (apply vc/get-vcf-retriever (cons ref-file vcf-files))]
          (vc/write-vcf-w-template (first vcf-files) {:out out-file}
                                   (map (comp (get-rep-vc vc-getter vcf-files (:names options)) parse-isec-line)
                                        (line-seq rdr))
                                   :header-update-fn (apply vc/merge-headers vcf-files)
                                   :new-md #{(VCFInfoHeaderLine. "CALLERS" VCFHeaderLineCount/UNBOUNDED
                                                                 VCFHeaderLineType/String
                                                                 "Individual caller support")}))
        (eprep/bgzip-index-vcf out-file)
                                        ;(fsp/remove-path work-dir)
        ))
    out-file))

(defn- usage [options-summary]
  (->> ["Ensemble calling for samples: combine multiple VCF caller outputs into a single callset."
        ""
        "Usage: bcbio-variation-recall ensemble [options] out-file ref-file [<vcf-files or list-files>]"
        ""
        "   out-file:   bgzipped VCF file to write merged output to"
        "   ref-file:   FASTA format genome reference file"
        "  <remaining>: VCF files to include for building a final ensemble callset."
        "               Specify on the command line or as text files containing paths to files."
        "               VCFs can be single or multi-sample."
        "               The input order of VCFs determines extraction preference in the final ensemble output."
        ""
        "Options:"
        options-summary]
       (string/join \newline)))

(defn -main [& args]
  (let [{:keys [options arguments errors summary]}
        (parse-opts args [["-c" "--cores CORES" "Number of cores to use" :default 1
                           :parse-fn #(Integer/parseInt %)]
                          ["-n" "--numpass NUMPASS" "Number of callers a variant should be present in to pass"
                           :default 2 :parse-fn #(Integer/parseInt %)]
                          [nil "--names NAMES"
                           "Comma separated list of names corresponding to VCFs for annotating output"
                           :default "" :parse-fn #(string/split % #",")]
                          [nil "--nofiltered" "Remove filtered variants before performing ensemble calls"]
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
