(ns bcbio.align.greads
  "Generic support for reads in multiple formats, transparently handles BAM/CRAM"
  (:require [bcbio.align.bam :as bam]
            [bcbio.align.cram :as cram]
            [bcbio.align.greads :as greads]
            [bcbio.run.itx :as itx]
            [bcbio.run.fsp :as fsp]
            [bcbio.variation.ensemble.prep :as eprep]
            [me.raynes.fs :as fs]))

(defmulti subset-in-region
  "Subset a read file to a region, handling BAM (no op) and CRAM"
  (fn [fname ref-file region tmp-dir]
    (keyword (subs (fs/extension fname) 1))))

(defmethod subset-in-region :bam
  ^{:doc "BAM format. No op: return original file"}
  [fname ref-file region tmp-dir]
  (bam/do-index fname)
  fname)

(defmethod subset-in-region :cram
  ^{:doc "CRAM subset, extract region of interest into temporary directory"}
  [fname ref-file region tmp-dir]
  (cram/do-index fname)
  (let [out-file (fsp/add-file-part fname (eprep/region->safestr region) tmp-dir ".bam")]
    (itx/run-cmd out-file
                 "scramble -I cram -O bam -r ~{ref-file} "
                 "-R ~{(eprep/region->samstr region)} ~{fname} ~{out-file}")
    (bam/do-index out-file)
    out-file))
