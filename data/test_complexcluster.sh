#!/bin/sh -e
#TODO: maybe change file name into filtercomplex.sh
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
        FNR==NR{chain_len[$1]=$3;next}
        {
            if (!($3 in off_arr)) {
                off_arr[$3]=OFFSET
            }
            cmpl_len[$3]+=chain_len[$1];OFFSET+=chain_len[$1]
        }
        END {
            for (cmpl in off_arr) {
                print cmpl"\t"off_arr[cmpl]"\t"cmpl_len[cmpl]
            }
        }' "${1}.index" "${1}.lookup" > "${2}.index"
    ln -s "$(abspath "${1}")" "${2}.0"
    cp "${1}.dbtype" "${2}.dbtype"
}

# check number of input variables
[ "$#" -ne 3 ] && echo "Please provide <queryDB> <targetDB> <clustDB> <tmpDir>" && exit 1;

# TODO : replace TMP_PATH
FILTALN="${QUERY}_filtcomp"
# DOING : filtercomplex
if notExists "${TMP_PATH}/${FILTALN}.dbtype"; then
    # shellcheck disable=SC2086
    $MMSEQS filtercomplex ${QUERY} ${TARGET} ${TMP_PATH}/${FILTALN}  ${FILTERCOMPLEX_PAR} \
        || fail "FilterComplex died"
fi

# FIXME : softlink source to complexDB
if notExists "${TMP_PATH}/cmpl_db.dbtype"; then
    buildCmplDb "${SOURCE}" "${TMP_PATH}/cmpl_db"
fi

INPUT="${TMP_PATH}/cmpl"
# FIXME : clust
if notExists "${TMP_PATH}/clu.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" clust "${INPUT}" ${TMP_PATH}/${FILTALN} "${RESULT}" ${CLUSTER_PAR} \
        || fail "Clustering died"
fi

# TODO : remove tmp
if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/${FILTALN}" ${VERBOSITY_PAR}
    "$MMSEQS" rmdb "${TMP_PATH}/cmpl_db" ${VERBOSITY_PAR}
    rm -rf ${TMP_PATH}/complexcluster.sh
