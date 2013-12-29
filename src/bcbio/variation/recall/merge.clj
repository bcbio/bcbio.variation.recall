(ns bcbio.variation.recall.merge
  "Merge multiple VCF files together, running in parallel over genomic regions."
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
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
  (fsp/add-file-part final-file (eprep/region->safestr region)
                     (fsp/safe-mkdir (io/file work-dir (get region :chrom "nochrom")))))

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
        vcf-file-str (string/join " " vcf-files)]
    (if (= 1 (count vcf-files))
      (io/copy (io/file (first vcf-files)) (io/file out-file))
      (itx/run-cmd out-file
                   "bcftools merge -o ~{(vcfutils/bcftools-out-type out-file)} "
                   "-r ~{(eprep/region->samstr region)} ~{vcf-file-str} "
                   "> ~{out-file}"))
    (eprep/bgzip-index-vcf out-file)))

(defmethod region-merge :bcftools
  ^{:doc "Merge VCFs within a region using bcftools."}
  [_ vcf-files region work-dir final-file]
  (let [group-size 200
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
          (doseq [ext ["" ".tbi"]]
            (.renameTo (io/file (str final-vcf ext)) (io/file (str out-file ext)))))))
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

(defmethod concatenate-vcfs :bgzip
  [vcf-files out-file config]
  (eprep/bgzip-index-vcf (concatenate-vcfs vcf-files (subs out-file 0 (.lastIndexOf out-file ".gz")) config)
                         :remove-orig? true))

(defmethod concatenate-vcfs :default
  [vcf-files out-file config]
  (let [str-vcf-files (string/join " " (rmap eprep/bgzip-index-vcf vcf-files (:cores config)))]
    (itx/run-cmd out-file
                 "vcfcat ~{str-vcf-files} > ~{out-file}")))

(defn combine-vcfs
  "Merge multiple VCF files together in parallel over genomic regions."
  [orig-vcf-files ref-file out-file config]
  (let [vcf-files (rmap eprep/bgzip-index-vcf orig-vcf-files (:cores config))
        _ (vcfutils/ensure-no-dup-samples vcf-files)
        merge-dir (fsp/safe-mkdir (io/file (fs/parent out-file) "merge"))
        region-bed (rsplit/group-pregions vcf-files ref-file merge-dir config)
        merge-parts (->> (rmap (fn [region]
                                 [(:i region) (region-merge :bcftools vcf-files region merge-dir out-file)])
                               (bed/reader region-bed)
                               (:cores config))
                         (map vec)
                         (sort-by first)
                         (map second))]
    (concatenate-vcfs merge-parts out-file config)))

(defn- is-vcf?
  [f]
  (with-open [rdr (vcfutils/pog-reader f)]
    (.startsWith (first (line-seq rdr)) "##fileformat=VCF")))

(defn- vcf-exists?
  [f]
  (or (and (fs/exists? f) (fs/file? f))
      (and (fs/exists? (str f ".gz")) (fs/file? (str f ".gz")))))

(defn- get-vcf-flex
  "Handle retrieving VCFs from either a single file or text-based list of files."
  [f]
  (cond
   (not (vcf-exists? f)) [f]
   (is-vcf? f) [f]
   :else (with-open [rdr (io/reader f)]
           (->> (line-seq rdr)
                (map string/trimr)
                (remove empty?)
                vec))))

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

(defn- error-msg [errors]
  (str "The following errors occurred while parsing your command:\n"
       (string/join \newline errors)))

(defn- exit [status msg]
  (println msg)
  (System/exit status))

(defn -main [& args]
  (let [{:keys [options arguments errors summary]}
        (parse-opts args [["-c" "--cores CORES" "Number of cores to use" :default 1
                           :parse-fn #(Integer/parseInt %)]
                          ["-h" "--help"]])]
    (cond
     (:help options) (exit 0 (usage summary))
     errors (exit 1 (error-msg errors))
     (= 0 (count arguments)) (exit 0 (usage summary))
     (< (count arguments) 3) (exit 1 (usage summary))
     :else (let [[out-file ref-file & vcf-inputs] arguments
                 vcf-files (mapcat get-vcf-flex vcf-inputs)
                 vcf-missing (remove vcf-exists? vcf-files)]
             (cond
              (not (empty? vcf-missing))
              (exit 1 (error-msg (cons "Input VCF files not found:" vcf-missing)))
              (or (not (fs/exists? ref-file)) (not (fs/file? ref-file)))
              (exit 1 (error-msg [(str "Reference file not found: " ref-file)]))
              :else
              (combine-vcfs vcf-files ref-file out-file options))))))
