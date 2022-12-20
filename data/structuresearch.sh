#!/bin/sh -e
# Assembler workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# 1. Finding exact $k$-mer matches.
if notExists "${TMP_PATH}/pref.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefilter "${QUERY_PREFILTER}" "${TARGET_PREFILTER}${INDEXEXT}" "${TMP_PATH}/pref" ${PREFILTER_PAR} \
        || fail "Kmer matching step died"
fi

# check if $ALIGNMENT_ALGO is tmalign
if [ "$ALIGNMENT_ALGO" = "tmalign" ]; then
    # 2. tm alignment
    if notExists "${TMP_PATH}/strualn.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" structurealign "${QUERY_ALIGNMENT}" "${TARGET_ALIGNMENT}${INDEXEXT}" "${TMP_PATH}/pref" "${TMP_PATH}/strualn" ${STRUCTUREALIGN_PAR} \
            || fail "Alignment step died"
    fi

    if notExists "${TMP_PATH}/aln.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" $ALIGNMENT_ALGO "${QUERY_ALIGNMENT}" "${TARGET_ALIGNMENT}${INDEXEXT}" "${TMP_PATH}/strualn" "${TMP_PATH}/aln" ${ALIGNMENT_PAR} \
            || fail "Alignment step died"
    fi

    if [ -n "$REMOVE_TMP" ]; then
        echo "Removing temporary files"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/strualn" ${VERBOSITY}
    fi
else
    # 2. Alignment
    if notExists "${TMP_PATH}/aln.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" $ALIGNMENT_ALGO "${QUERY_ALIGNMENT}" "${TARGET_ALIGNMENT}${INDEXEXT}" "${TMP_PATH}/pref" "${TMP_PATH}/aln" ${ALIGNMENT_PAR} \
            || fail "Alignment step died"
    fi
fi

# shellcheck disable=SC2086
"$MMSEQS" mvdb "${TMP_PATH}/aln" "${RESULTS}" ${VERBOSITY}

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/pref" ${VERBOSITY}
fi
