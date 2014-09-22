(ns bcbio.variation.test-erealign
  "Tests for ensemble variant consolidation with local realignment."
  (:require [bcbio.variation.ensemble.prep :as eprep]
            [bcbio.variation.ensemble.realign :as erealign]
            [bcbio.variation.ensemble.intersect :as eintersect]
            [bcbio.variation.recall.merge :as merge]
            [bcbio.variation.recall.split :as rsplit]
            [bcbio.variation.recall.square :as square]
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
               cram-file (str (io/file data-dir "NA12878-10.cram"))
               vcf-files (map #(str (io/file data-dir %))
                              ["NA12878-10-freebayes.vcf" "NA12878-10-gatk.vcf"
                               "NA12878-10-gatk-haplotype.vcf"])
               merge-vcf-file (str (io/file r-data-dir "NA12878-1-10-gatk-haplotype.vcf"))
               recall-bed (str (io/file r-data-dir "recall-regions.bed"))
               work-dir (str (io/file data-dir "work"))
               config {:cores 1 :caller :platypus}]
           (doseq [x (concat [work-dir]
                             (mapcat (fn [f]
                                       (map #(str (fsp/file-root f) %) [".vcf.gz" ".vcf.gz.tbi"
                                                                        "-passonly.vcf.gz" "-passonly.vcf.gz.tbi"]))
                                     (concat [merge-vcf-file] vcf-files)))]
             (fsp/remove-path x))
           ?form)))

(facts "Calculate ensemble set of variants from multiple inputs using realignment in a region."
  (let [region {:chrom "10" :start 250000 :end 399000}
        dirs {:ensemble work-dir}
        e-dir (fsp/safe-mkdir (io/file work-dir (eprep/region->safestr region)))
        out-file (str (io/file e-dir "recall-10_250000_399000.vcf.gz"))
        sample "NA12878-2"]
    (erealign/by-region sample region vcf-files bam-file ref-file dirs e-dir config) => out-file))

(facts "Calculate ensemble set of variants from multiple inputs using realignment over entire region."
  (let [out-file (str (io/file work-dir "NA12878-2-ensemble.vcf.gz"))]
    (erealign/ensemble-vcfs vcf-files [bam-file] ref-file out-file config) => out-file))

(facts "Calculate ensemble set of variants using intersection and n out of x counting."
  (let [out-file (str (io/file work-dir "NA12878-2-ensemble.vcf.gz"))
        econfig (assoc config :numpass 2)]
    (eintersect/ensemble-vcfs vcf-files ref-file out-file config) => out-file))

(facts "Identify split breakpoints for parallel execution"
  (let [out-file (str (io/file data-dir "work" "split" "NA12878-10-freebayes-combo-3-pregions.bed"))]
    (second (rsplit/group-pregions vcf-files ref-file work-dir config)) => out-file))

(facts "Merge multiple input files, running in parallel over small regions"
  (let [out-file (str (io/file data-dir "work" "NA12878-10-merge.vcf.gz"))]
    (merge/combine-vcfs [(first vcf-files) merge-vcf-file] ref-file out-file config) => out-file))

(facts "FreeBayes: Square off variant calls from multiple samples, creating merged final file."
  (let [out-file (str (io/file data-dir "work" "NA12878-10-square.vcf.gz"))
        fconfig (assoc config :caller :freebayes)]
    (square/combine-vcfs [(first vcf-files) merge-vcf-file] [bam-file bam-file]
                         ref-file out-file fconfig) => out-file))

(facts "Platypus: Square off variant calls from multiple samples, creating merged final file."
  (let [out-file (str (io/file data-dir "work" "NA12878-10-square.vcf.gz"))
        fconfig (assoc config :caller :platypus)]
    (square/combine-vcfs [(first vcf-files) merge-vcf-file] [bam-file bam-file]
                         ref-file out-file fconfig) => out-file))

(facts "samtools: Square off variant calls from multiple samples, creating merged final file."
  (let [out-file (str (io/file data-dir "work" "NA12878-10-square.vcf.gz"))
        fconfig (assoc config :caller :samtools)]
    (square/combine-vcfs [(first vcf-files) merge-vcf-file] [bam-file bam-file]
                         ref-file out-file fconfig) => out-file))

(facts "Squaring off with CRAM input files"
  (let [out-file (str (io/file data-dir "work" "NA12878-10-square.vcf.gz"))
        fconfig (-> config
                    (assoc :caller :freebayes)
                    ;(assoc :region recall-bed)
                    (assoc :region "10:300000-400000")
                    )]
    (square/combine-vcfs [(first vcf-files) merge-vcf-file] [cram-file cram-file]
                         ref-file out-file fconfig) => out-file))
