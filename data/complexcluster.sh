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
# Shift initial header DB into complex header DB
buildheadCmplDb() {
    awk -F"\t" '
    FNR==NR{
        split($2, parts, ".pdb")
        HEADS = parts[1]
        if (!($3 in cmplid)){
            cmplid[$3] = HEADS
            
        };
        next 
    }
    {
        split($1, part, ".pdb")
        COMP = substr(part[1], 2)
        if (!(COMP in head_arr)) {
            HEAD = substr(part[2], 4)
            head_arr[COMP] = HEAD
        }
    }
    END {
        for (cmpl in cmplid) {
            print cmplid[cmpl]".pdb", head_arr[cmplid[cmpl]]
        }
    }' "${1}.lookup" "${1}_h" > "${2}"
    cp "${1}.dbtype" "${2}.dbtype"
}

buildheadIndexCmplDb() {
    awk -F"\t" 'BEGIN {OFFSET=0}
        FNR==NR{
            if (!($3 in cmplchain)){
                cmplchain[$3]=$1
                };
            for (cmpl in cmplchain) {
                chaincmpl[cmplchain[cmpl]]=cmpl
                };
            next
            }
        {
            if (($1 in chaincmpl)) {
                off_arr[chaincmpl[$1]]=OFFSET
                OFFSET+=$3
                cmpl_len[chaincmpl[$1]]=$3
                }
        }
        END {
            for (cmpl in off_arr) {
                print cmpl"\t"off_arr[cmpl]"\t"cmpl_len[cmpl]
            }
        }' "${1}.lookup" "${1}_h.index" > "${2}.index"
}

# [ ! -d "$3" ] && echo "tmp directory $3 not found!" && mkdir -p "${TMP_PATH}";

if notExists "${TMP_PATH}/complex_result.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" complexsearch "${INPUT}" "${INPUT}" "${TMP_PATH}/complex_result" "${TMP_PATH}/complexsearch_tmp" ${COMPLEXSEARCH_PAR} \
        || fail "ComplexSearch died"
fi

if notExists "complex_filt.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" filtercomplex "${INPUT}" "${INPUT}" "${TMP_PATH}/complex_result" "${TMP_PATH}/complex_filt" ${FILTERCOMPLEX_PAR} \
        || fail "FilterComplex died"
fi

# shift query DB, .index, .dbtype
if notExists "${TMP_PATH}/complex_db.dbtype"; then    
    # build complex db as output
    buildCmplDb "${INPUT}" "${TMP_PATH}/complex_db"
fi
# Shift _h, _h.dbtype
if notExists "${TMP_PATH}/complex_db_h.dbtype"; then    
    # build complex header db as output
    buildheadCmplDb "${INPUT}" "${TMP_PATH}/complex_db_h"
fi

# Shift _h.index
if notExists "${TMP_PATH}/complex_db_h.index"; then   
    # build complex header.index as output 
    buildheadIndexCmplDb "${INPUT}" "${TMP_PATH}/complex_db_h"
fi

COMP="${TMP_PATH}/complex_db"

if notExists "${RESULT}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" clust "${COMP}" "${TMP_PATH}/complex_filt" "${RESULT}" ${CLUSTER_PAR} \
        || fail "Clustering died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/complex_filt" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/complex_result" ${VERBOSITY_PAR}
    rm -rf "${TMP_PATH}/complexsearch_tmp"
    rm -f "${TMP_PATH}/complexcluster.sh"
fi