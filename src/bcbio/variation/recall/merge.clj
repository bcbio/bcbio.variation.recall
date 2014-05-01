(ns bcbio.variation.recall.merge
  "Merge multiple VCF files together, running in parallel over genomic regions."
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.run.parallel :refer [rmap]]
            [bcbio.align.bed :as bed]
            [bcbio.variation.ensemble.prep :as eprep]
            [bcbio.variation.recall.clhelp :as clhelp]
            [bcbio.variation.recall.split :as rsplit]
            [bcbio.variation.recall.vcfutils :as vcfutils]
            [clojure.core.reducers :as r]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [clojure.tools.cli :refer [parse-opts]]
            [me.raynes.fs :as fs]))

(defn- region-merge-outfile
  "Build output file for regional merge inside of chromosome subdirectory."
  [region work-dir final-file]
  (fsp/add-file-part final-file (eprep/region->safestr region)
                     (fsp/safe-mkdir (io/file work-dir (get region :chrom "nochrom")))
                     ".vcf.gz"))

(defn- prep-vcf-region
  "Prepare a VCF file for merging, retrieving only the local region of interest."
  [vcf-file region tmp-dir]
  (let [out-file (fsp/add-file-part (string/replace vcf-file ".gz" "") (eprep/region->safestr region) tmp-dir)]
    (itx/run-cmd out-file
                 "tabix -h -p vcf ~{vcf-file}  ~{(eprep/region->samstr region)} "
                 "vcfcreatemulti > ~{out-file}")))

(defmulti region-merge
  "Perform merge of multiple sample VCFs within a given region."
  (fn [& args]
    (keyword (first args))))

(defn- run-bcftools
  "Run a bcftools merge on a (potential) subset of files in a temporary directory."
  [i vcf-files region out-dir base-file]
  (let [out-file (fsp/add-file-part base-file (format "%s-%s" (eprep/region->safestr region) i) out-dir)
        input-list (str (fsp/file-root out-file) "inputs.txt")]
    (if (= 1 (count vcf-files))
      (io/copy (io/file (first vcf-files)) (io/file out-file))
      (do
        (when (itx/needs-run? out-file)
          (spit input-list (string/join "\n" (map eprep/bgzip-index-vcf vcf-files))))
        (itx/run-cmd out-file
                     "bcftools merge -O ~{(vcfutils/bcftools-out-type out-file)} "
                     "-r ~{(eprep/region->samstr region)} `cat ~{input-list}` "
                     "> ~{out-file}")))
    (eprep/bgzip-index-vcf out-file)))

(defn move-vcf
  "Move a VCF file, also handling move of tabix index if it exists."
  [orig-file new-file]
  (doseq [ext ["" ".tbi"]]
    (when (fs/exists? (str orig-file ext))
      (.renameTo (io/file (str orig-file ext)) (io/file (str new-file ext))))))

(defmethod region-merge :bcftools
  ^{:doc "Merge VCFs within a region using bcftools."}
  [_ vcf-files region work-dir final-file]
  (let [group-size 5000
        out-file (region-merge-outfile region work-dir final-file)]
    (when (itx/needs-run? out-file)
      (itx/with-temp-dir [tmp-dir (fs/parent out-file)]
        (let [final-vcf (loop [work-files vcf-files
                               i 0]
                          (if (= 1 (count work-files))
                            (first work-files)
                            (let [merged (map-indexed (fn [j xs] (run-bcftools (+ i j) xs region tmp-dir out-file))
                                                      (partition-all group-size work-files))]
                              (recur merged (+ i (count merged))))))]
          (move-vcf final-vcf out-file))))
    out-file))

(defmethod region-merge :vcflib
  ^{:doc "Merge VCFs within a region using tabix and vcflib"}
  [_ vcf-files region work-dir final-file]
  (let [out-file (region-merge-outfile region work-dir final-file)]
    (when (itx/needs-run? out-file)
      (itx/with-temp-dir [tmp-dir (fs/parent out-file)]
        (let [prep-vcf-str (->> vcf-files
                                (map #(prep-vcf-region % region tmp-dir))
                                (string/join " "))
              bgzip-cmd (if (.endsWith out-file ".gz") "| bgzip -c" "")]
          (itx/run-cmd out-file
                       "vcfcombine ~{prep-vcf-str} | vcfcreatemulti ~{bgzip-cmd} > ~{out-file}"))))
    out-file))

(defmulti concatenate-vcfs
  "Concatenate VCF files in supplied order, handling bgzip and plain text."
  (fn [_ out-file _]
    (if (.endsWith out-file ".gz") :bgzip :default)))

(defn- do-concat
  [vcf-files out-file config]
  (let [input-list (str (fsp/file-root out-file) "-inputs.txt")
        bgzip-cmd (if (.endsWith out-file ".gz") "| bgzip -c" "")]
    (when (itx/needs-run? out-file)
      (spit input-list (string/join "\n" (rmap eprep/bgzip-index-vcf vcf-files (:cores config)))))
    (itx/run-cmd out-file
                 "vt concat `cat ~{input-list}` ~{bgzip-cmd} > ~{out-file}")))

(defmethod concatenate-vcfs :bgzip
  [vcf-files out-file config]
  (eprep/bgzip-index-vcf (do-concat vcf-files out-file config)))

(defmethod concatenate-vcfs :default
  [vcf-files out-file config]
  (do-concat vcf-files out-file config))

(defn prep-by-region
  "General functionality to split a set of VCFs into regions and apply a function, in parallel, to each."
  [f orig-vcf-files ref-file out-file config]
  (let [vcf-files (rmap #(eprep/bgzip-index-vcf % :remove-nopass? true) orig-vcf-files
                        (:cores config))
        merge-dir (fsp/safe-mkdir (io/file (fs/parent out-file) "merge"))
        region-bed (rsplit/group-pregions vcf-files ref-file merge-dir config)
        merge-parts (->> (rmap (fn [region]
                                 [(:i region) (f vcf-files region merge-dir)])
                               (bed/reader region-bed)
                               (:cores config))
                         (map vec)
                         (sort-by first)
                         (map second))]
    (concatenate-vcfs merge-parts out-file config)))

(defn combine-vcfs
  "Merge multiple VCF files together in parallel over genomic regions."
  [orig-vcf-files ref-file out-file config]
  (prep-by-region (fn [vcf-files region out-dir]
                    (vcfutils/ensure-no-dup-samples vcf-files)
                    (region-merge :bcftools vcf-files region out-dir out-file))
                  orig-vcf-files ref-file out-file config))

(defn- usage [options-summary]
  (->> ["Merge multiple VCF files together, running in parallel over genomic regions."
        ""
        "Usage: bcbio-variation-recall merge [options] out-file ref-file vcf-files"
        ""
        "  out-file:  VCF (or bgzipped VCF) file to write merged output to"
        "  ref-file:  FASTA format genome reference file"
        "  vcf-files: VCF files to merge. Can be specified on the command line "
        "             or as a text file containing paths to files for processing"
        ""
        "Options:"
        options-summary]
       (string/join \newline)))

(defn -main [& args]
  (let [{:keys [options arguments errors summary]}
        (parse-opts args [["-c" "--cores CORES" "Number of cores to use" :default 1
                           :parse-fn #(Integer/parseInt %)]
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
              (not (empty? (:bam arg-files)))
              (clhelp/exit 1 (clhelp/error-msg (cons "Did not expect BAM files:" (:bam arg-files))))
              (or (not (fs/exists? ref-file)) (not (fs/file? ref-file)))
              (clhelp/exit 1 (clhelp/error-msg [(str "Reference file not found: " ref-file)]))
              :else
              (combine-vcfs (:vcf arg-files) ref-file out-file options))))))
