(ns bcbio.variation.ensemble.realign
  "Realignment based Ensemble calling approaches using inputs from multiple callers.
   Uses tools from the Marth lab and Erik Garrison to realign and recall given a
   set of possible variants in a genomic region."
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.variation.ensemble.prep :as eprep]
            [bcbio.variation.recall.merge :as merge]
            [bcbio.variation.recall.square :as square]
            [bcbio.variation.recall.vcfutils :as vcfutils]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [me.raynes.fs :as fs]))

(defmulti realign-and-call
  (fn [& args]
    (keyword (get (last args) :caller :freebayes))))

;; TODO: Incorporate verifyBamID and --contamination-estimates
(defmethod realign-and-call :freebayes
  ^{:doc "Perform realignment with glia and calling with freebayes in a region of interest."}
  [region union-vcf bam-file ref-file work-dir config]
  (let [out-file (str (io/file work-dir (str "recall-" (eprep/region->safestr region) ".vcf")))]
    (itx/run-cmd out-file
                 "samtools view -bu ~{bam-file} ~{(eprep/region->samstr region)} | "
                 "glia -Rr -w 1000 -S 200 -Q 200 -G 4 -f ~{ref-file} -v ~{union-vcf} | "
                 "freebayes -f ~{ref-file} --variant-input ~{union-vcf} "
                 "--min-mapping-quality 1 --min-base-quality 3 --stdin | "
                 ;"vcffilter -f 'DP > 4' -f 'QUAL > 20.0' -t PASS -F QualDepthFilter | "
                 "vcfallelicprimitives -t MNP > ~{out-file}")
    out-file))

(defmethod realign-and-call :platypus
  ^{:doc "Perform realignment and recalling with platypus"}
  [region union-vcf bam-file ref-file work-dir config]
  (let [out-file (str (io/file work-dir (str "recall-" (eprep/region->safestr region) ".vcf.gz")))
        filters ["FR[*] <= 0.5 && TC < 4 && %QUAL < 20",
                 "FR[*] <= 0.5 && TC < 13 && %QUAL < 10",
                 "FR[*] > 0.5 && TC < 4 && %QUAL < 20"]
        filter_str (string/join " | " (map #(format "bcftools filter -e '%s' 2> /dev/null" %) filters))]
    (when (itx/needs-run? out-file)
      (itx/run-cmd out-file
                   "platypus callVariants --bamFiles=~{bam-file} --regions=~{(eprep/region->samstr region)} "
                   "--hapScoreThreshold 10 --scThreshold 0.99 --filteredReadsFrac 0.9 "
                   "--rmsmqThreshold 20 --qdThreshold 0 --abThreshold 0.00001 --minVarFreq 0.0 "
                   "--refFile=~{ref-file} --source=~{union-vcf} --assemble=1 "
                   "--logFileName /dev/null --verbosity=1 --output - "
                   "sed 's/\\tQ20\\t/\\tPASS\\t/' | "
                   "vt normalize -r ~{ref-file} -q - 2> /dev/null | vcfuniqalleles | "
                   "~{filter_str} | bgzip -c > ~{out-file}"))
    (eprep/bgzip-index-vcf out-file :remove-orig? true)))

(defn by-region
  "Realign and recall variants in a defined genomic region."
  [sample region vcf-files bam-file ref-file dirs e-dir config]
  (let [out-file (str (io/file e-dir (format "%s-%s.vcf.gz" sample (eprep/region->safestr region))))
        work-dir (str (.getParentFile (io/file out-file)))
        prep-dir (fsp/safe-mkdir (io/file work-dir (format "recall-%s-prep" (eprep/region->safestr region))))
        union-dir (fsp/safe-mkdir (io/file work-dir "union"))
        prep-inputs (map-indexed
                     (fn [i x]
                       (let [out-file (str (io/file prep-dir
                                                    (format "%s-%s-%s.vcf.gz"
                                                            sample (eprep/region->safestr region) i)))]
                         (square/subset-sample-region x sample region out-file)))
                         vcf-files)
        union-vcf (eprep/create-union :gatk prep-inputs ref-file region union-dir)]
    (realign-and-call region union-vcf bam-file ref-file work-dir config)))

(defn by-region-multi
  "Realign and recall variants in a region, handling multiple sample inputs."
  [vcf-files bam-files region ref-file dirs out-file config]
  (let [region-e-dir (fsp/safe-mkdir (io/file (:ensemble dirs) (get region :chrom "nochrom")
                                              (eprep/region->safestr region)))
        region-merge-dir (fsp/safe-mkdir (str region-e-dir "-merge"))
        final-vcfs (->> vcf-files
                        (mapcat #(for [s (vcfutils/get-samples %)] [s %]))
                        (group-by first)
                        (map (fn [[sample s-vcfs]]
                               (by-region sample region (map second s-vcfs) (get bam-files sample)
                                          ref-file dirs region-e-dir config))))]
    (merge/region-merge :gatk final-vcfs region ref-file region-merge-dir out-file)))

(defn ensemble-vcfs
  "Combine VCF files with squaring off by recalling at uncalled variant positions."
  [orig-vcf-files bam-files ref-file out-file config]
  (let [dirs {:ensemble (fsp/safe-mkdir (io/file (fs/parent out-file) "ensemble"))
              :inprep (fsp/safe-mkdir (io/file (fs/parent out-file) "inprep"))}
        bam-map (square/sample-to-bam-map bam-files ref-file (:inprep dirs))]
    (merge/prep-by-region (fn [vcf-files region merge-dir]
                            (by-region-multi vcf-files bam-map region ref-file
                                             (assoc dirs :merge merge-dir) out-file config))
                          orig-vcf-files ref-file out-file config)))
