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

(defn- region-merge
  "Perform merge of multiple sample VCFs within a given region."
  [vcf-files region work-dir final-file]
  (let [out-file (fsp/add-file-part final-file (eprep/region->safestr region)
                                    work-dir)
        vcf-file-str (string/join " " vcf-files)]
    (itx/run-cmd out-file
                 "bcftools merge -o ~{(vcfutils/bcftools-out-type out-file)} "
                 "-r ~{(eprep/region->samstr region)} ~{vcf-file-str} "
                 "> ~{out-file}")))

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
  [orig-vcf-files ref-file work-dir config]
  (vcfutils/ensure-no-dup-samples orig-vcf-files)
  (let [vcf-files (rmap eprep/bgzip-index-vcf orig-vcf-files (:cores config))
        region-bed (rsplit/group-pregions vcf-files ref-file work-dir)
        merge-dir (fsp/safe-mkdir (io/file (fs/parent region-bed) "merge"))
        out-file (fsp/add-file-part (fsp/remove-file-part region-bed "pregions")
                                    "merge" nil ".vcf.gz")
        merge-parts (->> (rmap (fn [region]
                                 [(:i region) (region-merge vcf-files region merge-dir out-file)])
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

(defn- get-vcf-flex
  "Handle retrieving VCFs from either a single file or text-based list of files."
  [f]
  (cond
   (or (not (fs/exists? f)) (not (fs/file? f))) [f]
   (is-vcf? f) [f]
   :else (with-open [rdr (io/reader f)]
           (->> (line-seq rdr)
                (map string/trimr)
                (remove empty?)
                vec))))

(defn- usage [options-summary]
  (->> ["Merge multiple VCF files together, running in parallel over genomic regions."
        ""
        "Usage: merge [options] out-file ref-file vcf-files-or-list"
        ""
        "vcf-files-or-list can be VCFs to merge specified on the command line or a text"
        "file containing the names of VCFs."
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
        (parse-opts args [["-c" "--cores" "Number of cores to use" :default 1
                           :parse-fn #(Integer/parseInt %)]
                          ["-h" "--help"]])]
    (cond
     (:help options) (exit 0 (usage summary))
     errors (exit 1 (error-msg errors))
     (= 0 (count arguments)) (exit 0 (usage summary))
     (< (count arguments) 3) (exit 1 (usage summary))
     :else (let [[out-file ref-file & vcf-inputs] arguments
                 vcf-files (mapcat get-vcf-flex vcf-inputs)
                 vcf-missing (remove #(and (fs/exists? %) (fs/file? %)) vcf-files)]
             (cond
              (not (empty? vcf-missing))
              (exit 1 (error-msg (cons "Input VCF files not found:" vcf-missing)))
              (or (not (fs/exists? ref-file)) (not (fs/file? ref-file)))
              (exit 1 (error-msg [(str "Reference file not found: " ref-file)]))
              :else
              (combine-vcfs vcf-files ref-file out-file options))))))
