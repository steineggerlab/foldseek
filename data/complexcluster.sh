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
[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[   -f "$2.dbtype" ] && echo "$2.dbtype already exists!" && exit 1;
[ ! -d "$3" ] && echo "tmp directory $3 not found!" && mkdir -p "$3";

INPUT="$1"
# DOING : complexsearch
if notExists "${TMP_PATH}/complex_result.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" complexsearch "${INPUT}" "${INPUT}" "${TMP_PATH}/complexsearch_aln" "${TMP_PATH}/complexsearch_tmp" ${COMPLEXSEARCH_PAR} \
        || fail "ComplexSearch died"
fi

# DOING : filtercomplex
if notExists "${RESULT}_filt.dbtype"; then
    # shellcheck disable=SC2086
    $MMSEQS filtercomplex "${INPUT}" "${INPUT}" "${TMP_PATH}/complexsearch_aln" "${RESULT}_filt" ${FILTERCOMPLEX_PAR} \
        || fail "FilterComplex died"
fi

# DOING : softlink source to complexDB
if notExists "${TMP_PATH}/cmpl_db.dbtype"; then
    buildCmplDb "${INPUT}" "${TMP_PATH}/cmpl_db"
fi

INPUTT="${TMP_PATH}/cmpl_db"

# DOING : clust
if notExists "${RESULT}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" clust "${INPUTT}" "${RESULT}_filt" "${RESULT}" ${CLUSTER_PAR} \
        || fail "Clustering died"
fi
 

# DOING: Remove tmp
if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${RESULT}_filt" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/complexsearch_aln" ${VERBOSITY_PAR}
    rm -rf "${TMP_PATH}/complexsearch_tmp"
    rm -f "${TMP_PATH}/easycomplexcluster.sh"
fi