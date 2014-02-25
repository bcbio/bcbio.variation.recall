(ns bcbio.variation.ensemble.realign
  "Realignment based Ensemble calling approaches using inputs from multiple callers.
   Uses tools from the Marth lab and Erik Garrison to realign and recall given a
   set of possible variants in a genomic region."
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.variation.ensemble.prep :as eprep]
            [bcbio.variation.recall.clhelp :as clhelp]
            [bcbio.variation.recall.merge :as merge]
            [bcbio.variation.recall.square :as square]
            [bcbio.variation.recall.vcfutils :as vcfutils]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [clojure.tools.cli :refer [parse-opts]]
            [me.raynes.fs :as fs]))

(defmulti realign-and-call
  (fn [& args]
    (keyword (get (last args) :caller :freebayes))))

;; TODO: Incorporate verifyBamID and --contamination-estimates
(defmethod realign-and-call :freebayes
  ^{:doc "Perform realignment with glia and calling with freebayes in a region of interest."}
  [region union-vcf bam-file ref-file work-dir config]
  (let [out-file (str (io/file work-dir (str "recall-" (eprep/region->safestr region) ".vcf")))]
    (itx/run-cmd out-file
                 "samtools view -bu ~{bam-file} ~{(eprep/region->samstr region)} | "
                 "glia -Rr -w 1000 -S 200 -Q 200 -G 4 -f ~{ref-file} -v ~{union-vcf} | "
                 "freebayes -f ~{ref-file} --variant-input ~{union-vcf} "
                 "--min-mapping-quality 1 --min-base-quality 3 --stdin | "
                 ;"vcffilter -f 'DP > 4' -f 'QUAL > 20.0' -t PASS -F QualDepthFilter | "
                 "vcfallelicprimitives -t MNP > ~{out-file}")
    out-file))

(defmethod realign-and-call :platypus
  ^{:doc "Perform realignment and recalling with platypus"}
  [region union-vcf bam-file ref-file work-dir config]
  (let [out-file (str (io/file work-dir (str "recall-" (eprep/region->safestr region) ".vcf.gz")))
        raw-out-file (string/replace out-file ".gz" "")]
    (when (itx/needs-run? out-file)
      (itx/run-cmd raw-out-file
                   "platypus callVariants --bamFiles=~{bam-file} --regions=~{(eprep/region->samstr region)} "
                   "--hapScoreThreshold 10 --scThreshold 0.99 "
                   "--refFile=~{ref-file} --source=~{union-vcf} --assemble=1 "
                   "--logFileName /dev/null --verbosity=1 --output ~{raw-out-file}")
      (square/platypus-filter raw-out-file))
    (eprep/bgzip-index-vcf raw-out-file :remove-orig? true)))

(defn by-region
  "Realign and recall variants in a defined genomic region."
  [sample region vcf-files bam-file ref-file dirs config]
  (let [out-file (str (io/file (fsp/safe-mkdir (io/file (:ensemble dirs) (get region :chrom "nochrom")))
                               (format "%s-%s.vcf.gz" sample (eprep/region->safestr region))))
        work-dir (str (.getParentFile (io/file out-file)))
        prep-dir (fsp/safe-mkdir (io/file work-dir (format "recall-%s-prep" (eprep/region->safestr region))))
        union-dir (fsp/safe-mkdir (io/file work-dir "union"))
        prep-inputs (map-indexed
                     (fn [i x]
                       (let [out-file (str (io/file prep-dir
                                                    (format "%s-%s-%s.vcf.gz"
                                                            sample (eprep/region->safestr region) i)))]
                         (square/subset-sample-region x sample region out-file)))
                         vcf-files)
        union-vcf (eprep/create-union prep-inputs ref-file region union-dir)]
    (realign-and-call region union-vcf bam-file ref-file work-dir config)))

(defn by-region-multi
  "Realign and recall variants in a region, handling multiple sample inputs."
  [vcf-files bam-files region ref-file dirs out-file config]
  (let [final-vcfs (->> vcf-files
                        (mapcat #(for [s (vcfutils/get-samples %)] [s %]))
                        (group-by first)
                        (map (fn [[sample s-vcfs]]
                               (by-region sample region (map second s-vcfs) (get bam-files sample)
                                          ref-file dirs config))))]
    (merge/region-merge :bcftools final-vcfs region (:merge dirs) out-file)))

(defn ensemble-vcfs
  "Combine VCF files with squaring off by recalling at uncalled variant positions."
  [orig-vcf-files bam-files ref-file out-file config]
  (let [dirs {:ensemble (fsp/safe-mkdir (io/file (fs/parent out-file) "ensemble"))}]
    (merge/prep-by-region (fn [vcf-files region merge-dir]
                            (by-region-multi vcf-files (square/sample-to-bam-map bam-files)
                                             region ref-file (assoc dirs :merge merge-dir) out-file config))
                          orig-vcf-files ref-file out-file config)))

(defn- usage [options-summary]
  (->> ["Ensemble calling for samples: combine multiple VCF caller outputs into a single callset."
        ""
        "Usage: bcbio-variation-recall ensemble [options] out-file ref-file [<vcf-files, bam-files, or list-files>]"
        ""
        "   out-file:   VCF (or bgzipped VCF) file to write merged output to"
        "   ref-file:   FASTA format genome reference file"
        "  <remaining>: VCF files to recall and BAM files for each sample. Can be specified "
        "               on the command line or as text files containing paths to files "
        "               for processing. VCFs can be single or multi-sample and BAMs can be in "
        "               any order but each VCF sample must have an associated BAM file."
        ""
        "Options:"
        options-summary]
       (string/join \newline)))

(defn -main [& args]
  (let [caller-opts #{:freebayes :platypus}
        {:keys [options arguments errors summary]}
        (parse-opts args [["-c" "--cores CORES" "Number of cores to use" :default 1
                           :parse-fn #(Integer/parseInt %)]
                          ["-m" "--caller CALLER" (str "Calling method to use: "
                                                       (string/join ", " (map name caller-opts)))
                           :default "platypus"
                           :parse-fn #(keyword %)
                           :validate [#(contains? caller-opts %)
                                      (str "Supported calling options: "
                                           (string/join ", " (map name caller-opts)))]]
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
              (ensemble-vcfs (:vcf arg-files) (:bam arg-files) ref-file out-file options))))))
