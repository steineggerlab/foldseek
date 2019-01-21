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
if notExists "${TMP_PATH}/pref"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" kmermatcher "$INPUT" "${TMP_PATH}/pref_$STEP" ${KMERMATCHER_PAR} \
        || fail "Kmer matching step died"
fi

# 2. tm alignment
if notExists "${TMP_PATH}/aln"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" align "$INPUT" "$INPUT" "${TMP_PATH}/pref" "${TMP_PATH}/aln" ${ALIGNMENT_PAR} \
        || fail "Ungapped alignment step died"
fi

  # 3. clust
if notExists "${TMP_PATH}/clust"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" align "$INPUT" "${TMP_PATH}/aln" "${TMP_PATH}/clu" ${CLUST_PAR} \
        || fail "Ungapped alignment step died"
fi

mv -f "${TMP_PATH}/clu" "$OUT_FILE" \
    || fail "Could not move result to $OUT_FILE"

mv -f "${TMP_PATH}/clu.index" "${OUT_FILE}.index" \
    || fail "Could not move result to ${OUT_FILE}.index"

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    rm -f "${TMP_PATH}/pref"
    rm -f "${TMP_PATH}/pref.index"
    rm -f "${TMP_PATH}/aln"
    rm -f "${TMP_PATH}/aln.index"
fi
