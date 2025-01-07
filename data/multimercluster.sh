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
    awk -F"\t" 'FNR==NR{
        split($2,parts,"_")
        multname=parts[1]
        for (j = 2; j < length(parts); j++) {
            if (j < length(parts)){
                multname=multname"_"
            }
            multname = multname parts[j]
        }
        multnameidx[multname]=$3; next
    } {
        split($1,words," ")
        split(words[1],part,"_")
        output_string=part[1]
        for (j = 2; j < length(part); j++) {
            if (j < length(part)){
                output_string=output_string"_"
            }
            output_string = output_string part[j]
        }
        if (output_string != current_key) {
            if (current_key != "") {
                print multnameidx[current_key] "\t" concatenated;
            }
            current_key = output_string;
            concatenated = $2;
        } else {
            concatenated = concatenated $2;
        }
    } END {
        if (current_key != "") {
            print multnameidx[current_key] "\t" concatenated;
        }
    }' "${2}.lookup" "${1}" > "${3}"
}

buldCmplhDb(){
    awk -F"\t" 'FNR==NR{
        split($2,parts,"_")
        multname=parts[1]
        for (j = 2; j < length(parts); j++) {
            if (j < length(parts)){
                multname=multname"_"
            }
            multname = multname parts[j]
        }
        multnameidx[multname]=$3; next
    }{
        split($2,words," ")
        split(words[1],part,"_")
        output_string=part[1]
        for (j = 2; j < length(part); j++) {
            if (j < length(part)){
                output_string=output_string"_"
            }
            output_string = output_string part[j]
        }
        headerstring=""
        for (k = 2; k < length(words)+1; k++) {
            headerstring = headerstring words[k]" "
        }
        print multnameidx[output_string]"\t"output_string" "headerstring
    }' "${2}.lookup" "${1}" | sort | uniq > "${3}"
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
    # build complex db as output
    # buildCmplDb "${INPUT}" "${TMP_PATH}/multimer_db"
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${INPUT}" "${INPUT}" "${TMP_PATH}/chain_db.tsv" --threads 1 ${VERBOSITY_PAR} \
        || fail "createtsv died"
    buildCmplDb "${TMP_PATH}/chain_db.tsv" "${INPUT}" "${TMP_PATH}/multimer.tsv"
    # shellcheck disable=SC2086
    "$MMSEQS" tsv2db "${TMP_PATH}/multimer.tsv" "${TMP_PATH}/multimer_db" --output-dbtype 0 ${VERBOSITY_PAR} \
        || fail "tsv2db died"
fi

BASEIN=$(basename "${INPUT}")
# Shift _h, _h.dbtype
if notExists "${TMP_PATH}/multimer_db_h.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" lndb "${INPUT}" "${TMP_PATH}/${BASEIN}tmp" ${VERBOSITY_PAR} \
        || fail "lndb died"

    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/${BASEIN}tmp_h" ${VERBOSITY_PAR} \
        || fail "rmdb died"

    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "${INPUT}.index" "${INPUT}_h" "${TMP_PATH}/${BASEIN}tmp_h" --subdb-mode 0 ${VERBOSITY_PAR} \
        || fail "createsubdb died"

    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/${BASEIN}tmp" "${TMP_PATH}/${BASEIN}tmp_h" "${TMP_PATH}/chain_db_h.tsv" --threads 1 ${VERBOSITY_PAR} \
        || fail "createtsv died"

    sort "${TMP_PATH}/chain_db_h.tsv" > "${TMP_PATH}/chain_db_h.tsvtmp"
    mv -f -- "${TMP_PATH}/chain_db_h.tsvtmp" "${TMP_PATH}/chain_db_h.tsv"
    buldCmplhDb "${TMP_PATH}/chain_db_h.tsv" "${INPUT}" "${TMP_PATH}/multimer_header.tsv"
    # shellcheck disable=SC2086
    "$MMSEQS" tsv2db "${TMP_PATH}/multimer_header.tsv" "${TMP_PATH}/multimer_db_h" --output-dbtype 12 ${VERBOSITY_PAR} \
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
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/${INPUT}tmp_h" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/${INPUT}tmp" ${VERBOSITY_PAR}
    rm -f -- "${TMP_PATH}/chain_db_h.tsv"
    rm -f -- "${TMP_PATH}/chain_db.tsv"
    rm -f -- "${TMP_PATH}/multimer.tsv"
    rm -f -- "${TMP_PATH}/multimer_header.tsv"
    rm -rf -- "${TMP_PATH}/multimersearch_tmp"
    rm -f -- "${TMP_PATH}/multimercluster.sh"
fi