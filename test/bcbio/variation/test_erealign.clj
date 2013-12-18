(ns bcbio.variation.test-erealign
  "Tests for ensemble variant consolidation with local realignment."
  (:require [bcbio.variation.ensemble.realign :as erealign]
            [bcbio.variation.recall.merge :as merge]
            [bcbio.variation.recall.split :as rsplit]
            [midje.sweet :refer :all]
            [clojure.java.io :as io]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]))

(background
 (around :facts
         (let [data-dir (str (io/file "." "test" "data" "ensemble"))
               r-data-dir (str (io/file "." "test" "data" "recall"))
               ref-file (str (io/file data-dir "chr10-start.fa"))
               bam-file (str (io/file data-dir "NA12878-10.bam"))
               vcf-files (map #(str (io/file data-dir %))
                              ["NA12878-10-freebayes.vcf" "NA12878-10-gatk.vcf"
                               "NA12878-10-gatk-haplotype.vcf"])
               merge-vcf-file (str (io/file r-data-dir "NA12878-1-10-gatk-haplotype.vcf"))
               work-dir (str (io/file data-dir "work"))]
           (doseq [x (concat [work-dir]
                             (map #(str % ".gz") vcf-files)
                             (map #(str % ".gz.tbi") vcf-files))]
             (fsp/remove-path x))
           ?form)))

(facts "Calculate ensemble set of variants from multiple inputs using realignment."
  (let [region {:chrom "10" :start 250000 :end 400000}
        out-file (str (io/file data-dir "work" "10_250000_400000-recall-10_250000_400000.vcf"))]
    (erealign/by-region region vcf-files bam-file ref-file work-dir) => out-file))

(facts "Identify split breakpoints for parallel execution"
  (let [out-file (str (io/file data-dir "work" "split" "NA12878-10-freebayes-combo-3-pregions.bed"))]
    (rsplit/group-pregions vcf-files ref-file work-dir) => out-file))

(facts "Merge multiple input files, running in parallel over small regions"
  (merge/combine-vcfs [(first vcf-files) merge-vcf-file] ref-file work-dir))
