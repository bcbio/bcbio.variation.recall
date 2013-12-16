(ns bcbio.variation.ensemble.realign
  "Realignment based Ensemble calling approaches using inputs from multiple callers.
   Uses tools from the Marth lab and Erik Garrison to realign and recall given a
   set of possible variants in a genomic region."
  (:require [bcbio.variation.ensemble.prep :as eprep]
            [bcbio.run.itx :as itx]
            [clojure.java.io :as io]
            [me.raynes.fs :as fs]))

(defn- realign-and-call
  "Perform realignment with glia and calling with freebayes in a region of interest.
   TODO: Incorporate verifyBamID and --contamination-estimates"
  [region union-vcf bam-file ref-file work-dir]
  (let [out-file (str (io/file work-dir (str "recall-" (eprep/region->safestr region) ".vcf")))]
    (itx/run-cmd out-file
                 "samtools view -bu ~{bam-file} ~{(eprep/region->samstr region)} | "
                 "glia -Rr -w 1000 -S 200 -Q 200 -G 4 -f ~{ref-file} -v ~{union-vcf} | "
                 "freebayes -f ~{ref-file} --haplotype-basis-alleles ~{union-vcf} "
                 "--min-mapping-quality 1 --min-base-quality 3 > ~{out-file}")
    out-file))

(defn by-region
  "Realign and recall variants in a defined genomic region."
  [region vcf-files bam-file ref-file work-dir]
  (let [region-dir (io/file work-dir (eprep/region->safestr region))
        _ (when-not (fs/exists? region-dir) (fs/mkdirs region-dir))
        prep-inputs (map #(eprep/norm-bgzip % region region-dir) vcf-files)
        union-vcf (eprep/create-union prep-inputs ref-file region region-dir)]
    (realign-and-call region union-vcf bam-file ref-file region-dir)))
