(ns bcbio.align.bed
  "Parse and manipulate tab-delimited BED files."
  (:require [clojure.core.protocols :as p]
            [clojure.core.reducers :as r]
            [clojure.string :as string]
            [clojure.java.io :as io]))

;; Provide a reducer-compatible list of regions in the BED file.
;; Hide resource management inside of the collection. Based on:
;; http://ce2144dc-f7c9-4f54-8fb6-7321a4c318db.s3.amazonaws.com/reducers.html#sec-2-2
(deftype BedReader [bed-file]
  p/CollReduce
  (p/coll-reduce [coll f]
    (p/coll-reduce coll f (f)))
  (p/coll-reduce [coll f init]
    (with-open [rdr (io/reader bed-file)]
      (loop [line-iter (line-seq rdr)
             v init]
        (if (empty? line-iter)
          v
          (let [r (f v (first line-iter))]
            (if (reduced? r)
              @r
              (recur (rest line-iter) r))))))))

(defn reader
  "High level bed reader that splits lines into chrom/start/end map"
  [bed-file]
  (r/map (fn [line]
           (let [[chrom start end] (take 3 (string/split line #"\t"))]
             {:chrom chrom :start (Integer/parseInt start) :end (Integer/parseInt end)}))
         (BedReader. bed-file)))
