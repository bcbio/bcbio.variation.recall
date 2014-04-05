(ns bcbio.align.cram
  "Manipulate CRAM files"
  (:require [bcbio.align.bam :as bam]
            [bcbio.align.gref :as gref]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [clojure.core.strint :refer [<<]]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [me.raynes.fs :as fs]))

(defn do-index
  "Perform indexing of a CRAM file inside a temporary directory"
  [cram-file]
  (let [cram-index (str cram-file ".crai")]
    (when (itx/needs-run? cram-index)
      (itx/with-tx-file [tx-cram-index cram-index]
        (let [tx-cram-file (fsp/file-root tx-cram-index)
              full-cram-file (str (fs/file cram-file))
              tmp-dir (str (fs/parent tx-cram-index))]
          (itx/check-run (<< "ln -s ~{full-cram-file} ~{tx-cram-file}"))
          (itx/check-run (<< "cram_index ~{tx-cram-file}")))))
    cram-index))

(defn sample-names
  "Identify sample names via conversion of CRAM to BAM to get the header"
  [cram-file ref-file]
  (do-index cram-file)
  (itx/with-temp-dir [tmp-dir (fs/parent cram-file)]
    (let [bam-file (fsp/add-file-part cram-file "header" tmp-dir ".bam")
          chrom (-> (gref/fasta-idx ref-file)
                    slurp
                    (string/split #"\n")
                    first
                    (string/split #"\t")
                    first)]
      (itx/run-cmd bam-file
                   "scramble -I cram -O bam -R ~{chrom}:1-1 ~{cram-file} ~{bam-file}")
      (bam/sample-names bam-file))))
