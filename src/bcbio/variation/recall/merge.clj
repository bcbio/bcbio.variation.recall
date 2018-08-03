(ns bcbio.variation.recall.merge
  "Merge multiple VCF files together, running in parallel over genomic regions."
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.run.clhelp :as clhelp]
            [bcbio.run.parallel :refer [rmap]]
            [bcbio.align.bed :as bed]
            [bcbio.variation.ensemble.prep :as eprep]
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
  (fsp/add-file-part final-file (eprep/region->safestr region) work-dir ".vcf.gz"))

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

(defn move-vcf
  "Move a VCF file, also handling move of tabix index if it exists."
  [orig-file new-file]
  (doseq [ext ["" ".tbi"]]
    (when (fs/exists? (str orig-file ext))
      (.renameTo (io/file (str orig-file ext)) (io/file (str new-file ext))))))

(defmethod region-merge :bcftools
  ^{:doc "Merge VCFs in a region using btools"}
  [_ vcf-files region ref-file final-file work-dir]
  (let [out-file (if work-dir (region-merge-outfile region work-dir final-file) final-file)
        input-list (str (fsp/file-root out-file) "-combineinputs.list")]
    (when (itx/needs-run? out-file)
      (spit input-list (string/join "\n" (map eprep/bgzip-index-vcf vcf-files))))
    (itx/with-temp-dir [tmp-dir (fs/parent out-file)]
      (itx/run-cmd out-file
                   "bcftools merge -O ~{(vcfutils/bcftools-out-type out-file)} "
                   "-r ~{(eprep/region->samstr region)} -l ~{input-list} "
                   "-o ~{out-file}"))
    (eprep/bgzip-index-vcf out-file)))

(defmethod region-merge :vcflib
  ^{:doc "Merge VCFs within a region using tabix and vcflib"}
  [_ vcf-files region ref-file work-dir final-file]
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
  (let [input-list (str (fsp/file-root out-file) "-concatinputs.list")
        bgzip-cmd (if (.endsWith out-file ".gz") "| bgzip -c" "")]
    (when (itx/needs-run? out-file)
      (spit input-list (string/join "\n" (rmap eprep/bgzip-index-vcf vcf-files (:cores config)))))
    (if (= 1 (count vcf-files))
      (itx/safe-copy (first vcf-files) out-file)
      (itx/run-cmd out-file
                   "bcftools concat --allow-overlaps --file-list ~{input-list} ~{bgzip-cmd} > ~{out-file}"))))

(defmethod concatenate-vcfs :bgzip
  [vcf-files out-file config]
  (eprep/bgzip-index-vcf (do-concat vcf-files out-file config)))

(defmethod concatenate-vcfs :default
  [vcf-files out-file config]
  (do-concat vcf-files out-file config))

(defn prep-by-region
  "General functionality to split a set of VCFs into regions and apply a function, in parallel, to each."
  [f orig-vcf-files ref-file out-file config]
  (let [in-dir (fsp/safe-mkdir (io/file (fs/parent out-file) "inprep"))
        bg-vcf-files (rmap #(eprep/bgzip-index-vcf % :remove-nopass? true :dir in-dir :orig-files orig-vcf-files)
                           orig-vcf-files (:cores config))
        merge-dir (fsp/safe-mkdir (io/file (fs/parent out-file) "merge"))
        [vcf-files region-bed] (rsplit/group-pregions bg-vcf-files ref-file merge-dir config)
        merge-parts (->> (rmap (fn [region]
                                 [(:i region) (f vcf-files region merge-dir)])
                               (bed/reader region-bed)
                               (:cores config))
                         (map vec)
                         (sort-by first)
                         (map second))]
    (when-not (empty? merge-parts)
      (concatenate-vcfs merge-parts out-file config))))

(defn combine-vcfs
  "Merge multiple VCF files together in parallel over genomic regions."
  [orig-vcf-files ref-file out-file config]
  (prep-by-region (fn [vcf-files region out-dir]
                    (vcfutils/ensure-no-dup-samples vcf-files)
                    (region-merge :bcftools vcf-files region ref-file out-file out-dir))
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
                          ["-r" "--region REGION"
                           "Genomic region to subset, in samtools format (chr1:100-200) or BED file"]
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
