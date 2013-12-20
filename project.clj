(defproject bcbio.variation.recall "0.0.1-SNAPSHOT"
  :description "Parallel merging, squaring off and ensemble calling for genomic variants."
  :url "https://github.com/chapmanb/bcbio.variation.recall"
  :license {:name "MIT" :url "http://www.opensource.org/licenses/mit-license.html"}
  :dependencies [[org.clojure/clojure "1.5.1"]
                 [org.clojure/tools.cli "0.3.0"]
                 [bcbio.run "0.0.1-SNAPSHOT"]]
  :plugins [[lein-midje "3.1.3"]]
  :profiles {:dev {:dependencies
                   [[midje "1.6.0"]]}
             :uberjar {:aot [bcbio.variation.recall.main]}}
  :main ^:skip-aot bcbio.variation.recall.main)
