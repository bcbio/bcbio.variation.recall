(ns bcbio.variation.recall.main
  (:require [bcbio.variation.recall.merge]
            [bcbio.variation.recall.square]
            [bcbio.variation.ensemble.intersect]
            [clojure.java.io :as io]
            [taoensso.timbre :as timbre])
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
            :doc (str "Perform squaring off for a set of called VCF files, recalling at no-call "
                      "positions in each sample.")}
   :ensemble {:main bcbio.variation.ensemble.intersect/-main
              :doc "Perform ensemble calling given multiple VCF callers for a single sample."}
   :version {:main version
             :doc "Print version"}})

(defn -main [& args]
  (if-let [to-run (get progs (keyword (first args)))]
    (do
      (try
        (apply (:main to-run) (rest args))
        (catch Exception e
          (timbre/error e)
          (shutdown-agents)
          (System/exit 1)))
      (shutdown-agents)
      (System/exit 0))
    (do
      (println "Parallel merging, squaring off and ensemble calling for genomic variants\n")
      (println "Commands:")
      (doseq [k (sort (keys progs))]
        (println (format "%-15s %s" (name k) (-> progs k :doc))))
      (System/exit 0))))
