# bcbio.variation.recall -- Create executable script

all:
		rm -f target/*.jar
		lein uberjar
	    cat bin/bcbio-variation-recall.template target/*-standalone.jar > bin/bcbio-variation-recall
	    chmod +x bin/bcbio-variation-recall
