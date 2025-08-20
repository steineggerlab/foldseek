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

getLookup() {
    awk 'FNR==NR{idx[$1]=1;next} {
        if ($1 in idx) {
            print $1"\t"$3
        }
    }' "${1}.index" "${1}.lookup" > "${2}.lookuptmp"
}

getSource() {
    awk 'FNR==NR{idx[$2]=1;next} {
        if ($1 in idx) {
            print
        }
    }' "${2}.lookuptmp" "${1}.source" > "${2}.sourcetmp"
}

buildIndex() {
    sort -nk2 "${1}.index" > "${1}.indextmp"
    awk 'FNR==NR{chainMult[$1]=$2; next} {
        print chainMult[$1]"\t"$2"\t"$3
    }' "${1}.lookuptmp" "${1}.indextmp" > "${1}.indextmp2"
     awk '!seen[$1]++ {
        sum[$1]=$3; off[$1]=$2; next
         } { sum[$1] += $3 } END{for(i in seen) {print i"\t"off[i]"\t"sum[i]
         }}' "${1}.indextmp2" > "${1}.index"
}

if notExists "${TMP_PATH}/multimer_result.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" multimersearch "${INPUT}" "${INPUT}" "${TMP_PATH}/multimer_result" "${TMP_PATH}/multimersearch_tmp" ${MULTIMERSEARCH_PAR} \
        || fail "multimerSearch died"
fi

if notExists "multimer_filt.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" filtermultimer "${INPUT}" "${INPUT}" "${TMP_PATH}/multimer_result" "${TMP_PATH}/multimer_filt" ${FILTERMULTIMER_PAR} \
        || fail "FilterMultimer died"
fi

# shift query DB, .index, .dbtype
if notExists "${TMP_PATH}/multimer_db.dbtype"; then
    getLookup "${INPUT}" "${TMP_PATH}/multimer_db"
    getSource "${INPUT}" "${TMP_PATH}/multimer_db"
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "${TMP_PATH}/multimer_db.lookuptmp" "${INPUT}" "${TMP_PATH}/multimer_db"  --subdb-mode 0 ${VERBOSITY_PAR} \
        || fail "createsubdb died"
    buildIndex "${TMP_PATH}/multimer_db"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/multimer_db_h" \
        || fail "rmdb died"
    # shellcheck disable=SC2086
    "$MMSEQS" tsv2db "${INPUT}.source" "${TMP_PATH}/multimer_db_h" --output-dbtype 12 ${VERBOSITY_PAR} \
        || fail "tsv2db died"
fi

COMP="${TMP_PATH}/multimer_db"

if notExists "${RESULT}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" clust "${COMP}" "${TMP_PATH}/multimer_filt" "${RESULT}" ${CLUSTER_PAR} \
        || fail "Clustering died"
    # shellcheck disable=SC2086
    "$MMSEQS" mvdb "${TMP_PATH}/multimer_filt_info" "${RESULT}_filt_info" ${VERBOSITY_PAR} \
        || fail "mv died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/multimer_filt" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/multimer_result" ${VERBOSITY_PAR}
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
