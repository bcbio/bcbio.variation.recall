(ns bcbio.variation.recall.main
  (:require [bcbio.variation.recall.merge]
            [bcbio.variation.recall.square]
            [clojure.java.io :as io])
  (:gen-class))

(defn version
  []
  (let [version (with-open [reader (-> "META-INF/maven/bcbio.variation.recall/bcbio.variation.recall/pom.properties"
                                       io/resource
                                       io/reader)]
                  (-> (doto (java.util.Properties.)
                        (.load reader))
                      (.getProperty "version")))]
    (println "bcbio.variation.recall" version)))

(def ^{:private true} progs
  {:merge {:main bcbio.variation.recall.merge/-main
           :doc "Merge multiple VCF files together, running in parallel over genomic regions."}
   :square {:main bcbio.variation.recall.square/-main
            :doc "Perform squaring off for a set of called VCF files, recalling at no-call positions in each sample."}
   :version {:main version
             :doc "Print version"}})

(defn -main [& args]
  (if-let [to-run (get progs (keyword (first args)))]
    (apply (:main to-run) (rest args))
    (do
      (println "Parallel merging, squaring off and ensemble calling for genomic variants\n")
      (println "Commands:")
      (doseq [k (sort (keys progs))]
        (println (format "%-15s %s" (name k) (-> progs k :doc)))))))
