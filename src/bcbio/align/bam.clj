(ns bcbio.align.bam
  "Manipulate BAM files, using the Picard samtools API"
  (:import [net.sf.samtools SAMFileReader SAMFileReader$ValidationStringency])
  (:require [clojure.java.io :as io]))

(defn sample-names
  "Retrieve samples represented in the BAM file."
  [bam-file]
  (SAMFileReader/setDefaultValidationStringency SAMFileReader$ValidationStringency/LENIENT)
  (with-open [in-bam (SAMFileReader. (io/file bam-file))]
    (set (map #(.getSample %) (-> in-bam .getFileHeader .getReadGroups)))))
