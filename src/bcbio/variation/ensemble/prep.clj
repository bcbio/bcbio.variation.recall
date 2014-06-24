(ns bcbio.variation.ensemble.prep
  "Prepare variant inputs for ensemble calling approaches.
   Creates normalized, bgzipped tabix indexed inputs from VCF files."
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [clojure.core.strint :refer [<<]]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [me.raynes.fs :as fs]))

(defn- tabix-index-vcf
  "Tabix index input VCF inside a transactional directory."
  [bgzip-file]
  (let [tabix-file (str bgzip-file ".tbi")]
    (when (or (itx/needs-run? tabix-file) (not (itx/up-to-date? tabix-file bgzip-file)))
      (itx/with-tx-file [tx-tabix-file tabix-file]
        (let [tx-bgzip-file (fsp/file-root tx-tabix-file)
              full-bgzip-file (str (fs/file bgzip-file))
              tmp-dir (str (fs/parent tx-bgzip-file))]
          (itx/check-run (<< "ln -s ~{full-bgzip-file} ~{tx-bgzip-file}"))
          (itx/check-run (<< "bcftools tabix -p vcf ~{tx-bgzip-file}")))))
    tabix-file))

(defn bgzip-index-vcf
  "Prepare a VCF file for positional query with bgzip and tabix indexing."
  [vcf-file & {:keys [remove-orig? remove-nopass? dir]}]
  (let [out-orig (str (fsp/file-root vcf-file)
                      (if remove-nopass? "-passonly" "")
                      ".vcf.gz")
        out-file (if dir (str (io/file dir (fs/base-name out-orig))) out-orig)]
    (if remove-nopass?
      (itx/run-cmd out-file "bcftools view -f 'PASS,.' ~{vcf-file} | bgzip -c > ~{out-file}")
      (itx/run-cmd out-file "bgzip -c ~{vcf-file} > ~{out-file}"))
    (when (and (not (.endsWith vcf-file ".gz")) remove-orig?)
      (fsp/remove-path vcf-file))
    (tabix-index-vcf out-file)
    out-file))

(defn region->samstr
  [region]
  (format "%s:%s-%s" (:chrom region) (inc (:start region)) (:end region)))

(defn region->safestr
  [region]
  (format "%s_%s_%s" (:chrom region) (:start region) (:end region)))

(defn region->freebayes
  [region]
  (format "%s:%s..%s" (:chrom region) (:start region) (:end region)))

(defn norm-bgzip
  "Normalize and bgzip/tabix index a VCF input file in a defined region."
  [vcf-file region out-dir]
  (let [prep-vcf-file (bgzip-index-vcf vcf-file)
        out-file (str (io/file out-dir (str (fs/base-name vcf-file) ".gz")))]
    (itx/run-cmd out-file
             "bcftools view -r ~{(region->samstr region)} ~{prep-vcf-file} | "
             "vcfallelicprimitives | "
             "bgzip -c /dev/stdin > ~{out-file}")
    (tabix-index-vcf out-file)
    out-file))

(defmulti create-union
  "Create a minimal union file with inputs from multiple variant callers in the given region."
  (fn [& args]
    (keyword (first args))))

(defmethod create-union :gatk
  ^{:doc "Use GATK CombineVariants to merge multiple input files with sites-only output"}
  [_ vcf-files ref-file region out-dir]
  (let [out-file (str (io/file out-dir (str "union-" (region->safestr region) ".vcf.gz")))
        variant-str (string/join " " (map #(str "--variant " (bgzip-index-vcf %)) vcf-files))]
    (itx/run-cmd out-file
                 "gatk-framework -Xms250m -Xmx2g -T CombineVariants -R ~{ref-file} "
                 "-L ~{(region->samstr region)} --out ~{out-file} "
                 "--minimalVCF --sites_only "
                 "--suppressCommandLineHeader --setKey null "
                 "-U LENIENT_VCF_PROCESSING --logging_level ERROR "
                 "~{variant-str}")
    (bgzip-index-vcf out-file)))

(defmethod create-union :bcftools
  ^{:doc "Use bcftools isec and custom awk command to handle merge of multiple files"}
  [_ vcf-files ref-file region out-dir]
  (let [out-file (str (io/file out-dir (str "union-" (region->safestr region) ".vcf.gz")))
        vcf-files-str (string/join " " vcf-files)
        vcf-header "echo -e '##fileformat=VCFv4.1\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO'"
        isec-to-vcf "awk -F'\\t' '{ print $1 FS $2 FS \".\" FS $3 FS $4 FS \".\" FS \".\" FS \".\"}'"]
    (itx/run-cmd out-file
                 "cat <(~{vcf-header}) "
                 "<(bcftools isec -n +1 -r ~{(region->samstr region)} ~{vcf-files-str} | ~{isec-to-vcf)) | "
                 "vcfcreatemulti | bgzip -c > ~{out-file}")
    (bgzip-index-vcf out-file :remove-orig? true)))
