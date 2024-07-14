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

buldCmplhDb(){
    awk -F"\t" 'BEGIN {INDEXVAL=0}
        {
            split($2,words," ")
            split(words[1],parts,"_")
            output_string=""
            for (j = 1; j < length(parts); j++) {
                output_string = output_string parts[j]
                if (j < length(parts)-1){
                    output_string=output_string"_" 
                }
            }
            headerstring=""
            for (k = 2; k < length(words)+1; k++) {
                headerstring = headerstring words[k]" "
            }
            if (!(output_string in gogo)){
                print INDEXVAL"\t"output_string" "headerstring
                INDEXVAL++
            }
            gogo[output_string]=1
            
        }' "${1}" > "${2}"
}


# [ ! -d "$3" ] && echo "tmp directory $3 not found!" && mkdir -p "${TMP_PATH}";

if notExists "${TMP_PATH}/complex_result.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" multimersearch "${INPUT}" "${INPUT}" "${TMP_PATH}/complex_result" "${TMP_PATH}/multimersearch_tmp" ${COMPLEXSEARCH_PAR} \
        || fail "multimerSearch died"
fi

if notExists "complex_filt.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" filtercomplex "${INPUT}" "${INPUT}" "${TMP_PATH}/complex_result" "${TMP_PATH}/complex_filt" "${TMP_PATH}/filtcov.tsv" ${FILTERCOMPLEX_PAR} \
        || fail "FilterComplex died"
fi

# shift query DB, .index, .dbtype
if notExists "${TMP_PATH}/complex_db.dbtype"; then    
    # build complex db as output
    buildCmplDb "${INPUT}" "${TMP_PATH}/complex_db"
fi

# Shift _h, _h.dbtype
if notExists "${TMP_PATH}/complex_db_h.dbtype"; then
    # # shellcheck disable=SC2086
    # "$MMSEQS" tsv2db "${INPUT}.source" "${TMP_PATH}/complex_db_header_tmp" ${VERBOSITY_PAR} \
    #     || fail "tsv2db died"
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${INPUT}" "${INPUT}_h" "${TMP_PATH}/chain_db_h.tsv" ${VERBOSITY_PAR} \
        || fail "createtsv died"
    buldCmplhDb "${TMP_PATH}/chain_db_h.tsv" "${TMP_PATH}/complex_header.tsv"
    # shellcheck disable=SC2086
    "$MMSEQS" tsv2db "${TMP_PATH}/complex_header.tsv" "${TMP_PATH}/complex_db_h" ${VERBOSITY_PAR} \
        || fail "tsv2db died"
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
    rm "${TMP_PATH}/complex_header.tsv"
    rm "${TMP_PATH}/complex_header.tsv_redundant"
    rm -rf "${TMP_PATH}/complexsearch_tmp"
    rm -f "${TMP_PATH}/complexcluster.sh"
fi