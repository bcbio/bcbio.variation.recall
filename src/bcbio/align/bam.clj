(ns bcbio.align.bam
  "Manipulate BAM files, using the Picard samtools API"
  (:import [htsjdk.samtools SAMFileReader ValidationStringency])
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [clojure.core.strint :refer [<<]]
            [clojure.java.io :as io]
            [me.raynes.fs :as fs]))

(defn do-index
  "Index BAM file"
  [bam-file]
  (let [out-file (str bam-file ".bai")]
    (when (or (itx/needs-run? out-file) (not (itx/up-to-date? out-file bam-file)))
      (itx/with-tx-file [tx-out-file out-file]
        (let [tx-bam-file (fsp/file-root tx-out-file)
              full-bam-file (str (fs/file bam-file))]
          (itx/check-run (<< "ln -s ~{full-bam-file} ~{tx-bam-file}"))
          (itx/check-run (<< "sambamba index ~{tx-bam-file}")))))))

(defn sample-names
  "Retrieve samples represented in the BAM file."
  [bam-file]
  (SAMFileReader/setDefaultValidationStringency ValidationStringency/LENIENT)
  (with-open [in-bam (SAMFileReader. (io/file bam-file))]
    (set (map #(.getSample %) (-> in-bam .getFileHeader .getReadGroups)))))
