(ns bcbio.variation.recall.split
  "Split variation calls into small genomic regions for parallel processing."
  (:require [bcbio.align.gref :as gref]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.variation.recall.vcfutils :as vcfutils]
            [clojure.java.io :as io]
            [schema.core :as s]))

(s/defn ^:always-validate vcf-breakpoints
  "Prepare BED file of non-variant regions in the input VCF as parallel breakpoints."
  ([vcf-file sample :- s/String ref-file work-dir]
     (let [merge-size 500
           split-dir (fsp/safe-mkdir (io/file work-dir sample "split"))
           sample-vcf (vcfutils/subset-to-sample vcf-file sample split-dir)
           out-file (fsp/add-file-part sample-vcf "splitpoints" split-dir ".bed")
           fai-file (gref/fasta-idx ref-file)]
       (itx/run-cmd out-file
                    "bedtools genomecov -i ~{sample-vcf} -g ~{fai-file} -bg | "
                    "bedtools merge -d ~{merge-size} > ~{out-file}")
       out-file))
  ([vcf-file ref-file work-dir]
     (vcf-breakpoints vcf-file (vcfutils/get-sample vcf-file) ref-file work-dir)))
