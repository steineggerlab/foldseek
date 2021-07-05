#!/bin/sh -e
# Assembler workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# check input variables
[ ! -n "${OUT_FILE}" ] && echo "Please provide OUT_FILE" && exit 1
[ ! -n "${TMP_PATH}" ] && echo "Please provide TMP_PATH" && exit 1

# check if files exists
[   -f "${OUT_FILE}" ] &&  echo "${OUT_FILE} exists already!" && exit 1
[ ! -d "${TMP_PATH}" ] &&  echo "tmp directory ${TMP_PATH} not found!" && mkdir -p "${TMP_PATH}"

# 1. Finding exact $k$-mer matches.
if notExists "${TMP_PATH}/pref.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefilter "$INPUT" "${TMP_PATH}/pref_$STEP" ${PREFILTER_PAR} \
        || fail "Kmer matching step died"
fi

# 2. tm alignment
if notExists "${TMP_PATH}/aln.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" align "$INPUT" "$INPUT" "${TMP_PATH}/pref" "${TMP_PATH}/aln" ${ALIGNMENT_PAR} \
        || fail "Ungapped alignment step died"
fi

"$MMSEQS" mvdb "${TMP_PATH}/aln" "$OUT_FILE"

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    rm -f "${TMP_PATH}/pref"
    rm -f "${TMP_PATH}/pref.index"
    rm -f "${TMP_PATH}/aln"
    rm -f "${TMP_PATH}/aln.index"
fi
