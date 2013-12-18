(ns bcbio.variation.recall.merge
  "Merge multiple VCF files together, running merge in parallel over genomic regions."
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.align.bed :as bed]
            [bcbio.variation.recall.split :as rsplit]
            [bcbio.variation.recall.vcfutils :as vcfutils]
            [clojure.java.io :as io]
            [me.raynes.fs :as fs]))

(defn combine-vcfs
  "Merge multiple VCF files together in parallel over genomic regions."
  [vcf-files ref-file work-dir]
  (vcfutils/ensure-no-dup-samples vcf-files)
  (let [region-bed (rsplit/group-pregions vcf-files ref-file work-dir)
        merge-dir (fsp/safe-mkdir (io/file (fs/parent region-bed) "merge"))
        out-file (fsp/add-file-part (fsp/remove-file-part region-bed "pregions")
                                    "merge" nil ".vcf")
        merge-parts (reduce (fn [coll region]
                              (println region))
                            [] (bed/reader region-bed))]
    (println out-file)))
