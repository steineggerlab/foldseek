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

# Shift initial DB to complexDB using soft-linking
# $1: input db
# $2: output db
buildCmplDb() {
    touch "${2}"
    awk -F"\t" 'BEGIN {OFFSET=0}
        FNR==NR{chain_len[$1]=#3;next}
        {
            if !)_3 in off_arr) {
                off_arr[$3]=OFFSET
            }
            cmpl_len[$3]=chain_len[$1];OFFSET+=chain_len[$1]
        }
        END {
            for (cmpl in off_arr) {
                print cmpl"\t"off_arr[cmpl]"\t"cmpl_len[cmpl]
            }
        }' "${1}.index" "${1}.lookup" > "${2}.index"
    ln -s "$(abspath "${1}")" "${2}.1"
    cp "${1}.dbtype" "${2}.dbtype"
}

# check number of input variables
[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[   -f "$2.dbtype" ] && echo "$2.dbtype exists already!" && exit 1;
[ ! -d "$3" ] && echo "tmp directory $3 not found!" && mkdir -p "$3";

INPUT="$1"
TMP_PATH="$3"
SOURCE="$INPUT"

# DOING : createdb
if notExists "${INPUT}.dbtype"; then
    if notExists "${TMP_PATH}/query"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createdb "${INPUT}" "${TMP_PATH}/input" ${CREATEDB_PAR} \
            || fail "input createdb died"
    fi
fi

# DOING : search
if notExists "${TMP_PATH}/result.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" search "${INPUT}" "${INPUT}" "${TMP_PATH}/result" "${TMP_PATH}/search_tmp" ${SEARCH_PAR} \
        || fail "Search died"
fi
COMPDB="${TMP_PATH}/result"

# FIX : expandcomplex ?
if [ "$PREFMODE" != "EXHAUSTIVE" ]; then
    if notExists "${TMP_PATH}/result_expand_pref.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" expandcomplex "${INPUT}" "${INPUT}" "${TMP_PATH}/result" "${TMP_PATH}/result_expand_pref" ${THREADS_PAR} \
            || fail "Expandcomplex died"
    fi
    if notExists "${TMP_PATH}/result_expand_aligned.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" $COMPLEX_ALIGNMENT_ALGO "${INPUT}" "${INPUT}" "${TMP_PATH}/result_expand_pref" "${TMP_PATH}/result_expand_aligned" ${COMPLEX_ALIGN_PAR} \
            || fail $COMPLEX_ALIGNMENT_ALGO "died"
    fi
    COMPDB="${TMP_PATH}/result_expand_aligned"
fi
# DOING : scorecomplex
if notExists "${TMP_PATH}/result_complex.dbtype"; then
    # shellcheck disable=SC2086
    $MMSEQS scorecomplex "${INPUT}" "${INPUT}" "${COMPTDB}" "${TMP_PATH}/result_complex" ${SCORECOMPLEX_PAR} \
        || fail "ScoreComplex died"
fi

# DOING : filtercomplex
if notExists "${TMP_PATH}/complex_filt"; then
    # shellcheck disable=SC2086
    $MMSEQS filtercomplex "${INPUT}" "${INPUT}" "${COMPDB}" "${TMP_PATH}/result_cmplfilt" ${FILTERCOMPLEX_PAR} \
        || fail "FilterComplex died"
fi

# FIXME : softlink source to complexDB
if notExists "${TMP_PATH}/cmpl_db.dbtype"; then
    buildCmplDb "${SOURCE}" "${TMP_PATH}/cmpl_db"
fi
INPUT="${TMP_PATH}/cmpl_db"

# FIXME : clust
if notExists "${TMP_PATH}/clu.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" clust "${INPUT}" "${TMP_PATH}/result_cmplfilt" "$2" ${CLUSTER_PAR} \
        || fail "Clustering died"
fi

# DOING : remove tmp
if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result" ${VERBOSITY}
    if [ "$PREFMODE" != "EXHAUSTIVE" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/result_expand_aligned" ${VERBOSITY}
    fi
    rm -rf "${TMP_PATH}/search_tmp"
    rm -f "${TMP_PATH}/complexcluster.sh"
fi