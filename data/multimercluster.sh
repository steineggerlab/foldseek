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

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [ -z "${1##*/*}" ]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

if notExists "${TMP_PATH}/multimer_result.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" multimersearch "${INPUT}" "${INPUT}" "${TMP_PATH}/multimer_result" "${TMP_PATH}/multimersearch_tmp" ${MULTIMERSEARCH_PAR} \
        || fail "multimerSearch died"
fi

COMP="${TMP_PATH}/multimer_result_query_multimerdb"

if notExists "${RESULT}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" clust "${COMP}" "${TMP_PATH}/multimer_result" "${RESULT}" ${CLUSTER_PAR} \
        || fail "Clustering died"
    # shellcheck disable=SC2086
    "$MMSEQS" setextendeddbtype "${RESULT}" --extended-dbtype 16 ${VERBOSITY_PAR} \
        || fail "setextendeddbtype died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/multimer_result" ${VERBOSITY_PAR} \
        || fail "rmdb died"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/multimer_result_query_multimerdb" ${VERBOSITY_PAR} \
        || fail "rmdb died"
    rm -f -- "${TMP_PATH}/chain_db_h.tsv"
    rm -f -- "${TMP_PATH}/chain_db.tsv"
    rm -f -- "${TMP_PATH}/multimer.tsv"
    rm -f -- "${TMP_PATH}/multimer_header.tsv"
    rm -rf -- "${TMP_PATH}/multimersearch_tmp"
    rm -f -- "${TMP_PATH}/multimer_db.lookuptmp"
    rm -f -- "${TMP_PATH}/multimer_db.sourcetmp"
    rm -f -- "${TMP_PATH}/multimer_db.indextmp2"
    rm -f -- "${TMP_PATH}/multimer_db.indextmp"
    rm -f -- "${TMP_PATH}/multimercluster.sh"
fi
