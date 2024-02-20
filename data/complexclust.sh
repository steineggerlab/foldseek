#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

exists() {
	[ -f "$1" ]
}

# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <sequenceDB> <clusterDB> <clusterDB> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[   -f "$3.dbtype" ] && echo "$3.dbtype exists already!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";

INPUTSEQ="$1"
INPUTCLUST=="$2" 
OUTCLUST="$3"
TMP_PATH="$4"

COMPLEXINPUT="${TMP_PATH}/${INPUTSEQ}_com"


"$MMSEQS" clust "$COMPLEXINPUT" "$INPUTCLUST" "$OUTCLUST" ${CLUSTER_PAR}

