(ns bcbio.align.bam
  "Manipulate BAM files, using the Picard samtools API"
  (:import [net.sf.samtools SAMFileReader SAMFileReader$ValidationStringency])
  (:require [bcbio.run.itx :as itx]
            [clojure.java.io :as io]))

(defn do-index
  "Index BAM file"
  [bam-file]
  (let [out-file (str bam-file ".bai")]
    (itx/run-cmd out-file
                 "sambamba index ~{bam-file} ~{out-file}")))

(defn sample-names
  "Retrieve samples represented in the BAM file."
  [bam-file]
  (SAMFileReader/setDefaultValidationStringency SAMFileReader$ValidationStringency/LENIENT)
  (with-open [in-bam (SAMFileReader. (io/file bam-file))]
    (set (map #(.getSample %) (-> in-bam .getFileHeader .getReadGroups)))))
