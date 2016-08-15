(defproject bcbio.variation.recall "0.1.7"
  :description "Parallel merging, squaring off and ensemble calling for genomic variants."
  :url "https://github.com/chapmanb/bcbio.variation.recall"
  :license {:name "MIT" :url "http://www.opensource.org/licenses/mit-license.html"}
  :dependencies [[org.clojure/clojure "1.7.0"]
                 [org.clojure/tools.cli "0.3.0"]
                 [ordered "1.3.2"]
                 [version-clj "0.1.0"]
                 [de.kotka/lazymap "3.1.1"]
                 [bcbio.run "0.0.5"]
                 [com.github.broadinstitute/picard "1.140"]
                 [com.github.samtools/htsjdk "1.140"]]
  :plugins [[lein-midje "3.1.3"]]
  :profiles {:dev {:dependencies
                   [[midje "1.6.3"]]}
             :uberjar {:aot [bcbio.variation.recall.main]}}
  :main ^:skip-aot bcbio.variation.recall.main)
