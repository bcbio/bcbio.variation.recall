(ns bcbio.align.gref
  "Handle genomic references in FASTA format."
  (:require [bcbio.run.itx :as itx]))

(defn fasta-idx
  "Create idx file for fasta reference input."
  [fasta-file]
  (let [out-file (str fasta-file ".fai")]
    (itx/run-cmd out-file
                 "samtools faidx ~{fasta-file}")))
