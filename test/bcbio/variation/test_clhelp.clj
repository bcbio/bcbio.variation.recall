(ns bcbio.variation.test-clhelp
  "Test support functionality for writing useful command lines."
  (:require [bcbio.variation.recall.clhelp :as clhelp]
            [midje.sweet :refer :all]
            [clojure.java.io :as io]))

(facts "Extract command lines containing mixture of VCF, BAM and list files."
  (let [data-dir (str (io/file "." "test" "data" "ensemble"))
        r-data-dir (str (io/file "." "test" "data" "recall"))
        real-vcfs (map #(str (io/file data-dir %))
                       ["NA12878-10-freebayes.vcf" "NA12878-10-gatk.vcf"
                        "NA12878-10-gatk-haplotype.vcf"])
        real-bams [(str (io/file data-dir "NA12878-10.bam"))]
        wrong-files ["not-a-file.vcf" "not-real.bam"]]
    (clhelp/vcf-bam-args (concat real-vcfs real-bams)) => {:vcf real-vcfs :bam real-bams}
    (clhelp/vcf-bam-args (concat real-vcfs real-bams wrong-files)) => {:vcf real-vcfs :bam real-bams
                                                                       :missing wrong-files}
    (clhelp/vcf-bam-args [(str (io/file r-data-dir "arg-list-good.txt"))]) => {:vcf real-vcfs :bam real-bams}
    (clhelp/vcf-bam-args [(str (io/file r-data-dir "arg-list-bad.txt"))]) => {:vcf real-vcfs :bam real-bams
                                                                              :missing wrong-files}))
