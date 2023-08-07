#!/bin/sh -e
# Assembler workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
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

fake_pref() {
    QDB="$1"
    TDB="$2"
    RES="$3"
    # create link to data file which contains a list of all targets that should be aligned
    ln -s "$(abspath "${TDB}.index")" "${RES}"
    # create new index repeatedly pointing to same entry
    INDEX_SIZE="$(wc -c < "${TDB}.index")"
    awk -v size="$INDEX_SIZE" '{ print $1"\t0\t"size; }' "${QDB}.index" > "${RES}.index"
    # create dbtype (7)
    awk 'BEGIN { printf("%c%c%c%c",7,0,0,0); exit; }' > "${RES}.dbtype"
}

# 1. Finding exact $k$-mer matches.
if notExists "${TMP_PATH}/pref.dbtype"; then
    if [ "$PREFMODE" = "EXHAUSTIVE" ]; then
        fake_pref "${QUERY_PREFILTER}" "${TARGET_PREFILTER}" "${TMP_PATH}/pref"
    elif [ "$PREFMODE" = "UNGAPPED" ]; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" ungappedprefilter "${QUERY_PREFILTER}" "${TARGET_PREFILTER}${INDEXEXT}" "${TMP_PATH}/pref" ${UNGAPPEDPREFILTER_PAR} \
            || fail "Ungapped prefilter matching step died"
    else
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" prefilter "${QUERY_PREFILTER}" "${TARGET_PREFILTER}${INDEXEXT}" "${TMP_PATH}/pref" ${PREFILTER_PAR} \
            || fail "Kmer matching step died"
    fi
fi

# check if $ALIGNMENT_ALGO is tmalign
if [ "$ALIGNMENT_ALGO" = "tmalign" ]; then
    # 2. tm alignment
    INTERMEDIATE="${TMP_PATH}/strualn"
    if notExists "${TMP_PATH}/strualn.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" structurealign "${QUERY_ALIGNMENT}" "${TARGET_ALIGNMENT}${INDEXEXT}" "${TMP_PATH}/pref" "${INTERMEDIATE}" ${STRUCTUREALIGN_PAR} \
            || fail "Alignment step died"
    fi

    if [ -n "${EXPAND}" ]; then
        if notExists "${TMP_PATH}/strualn_expanded.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" mergeresultsbyset  "${INTERMEDIATE}" "${TARGET_ALIGNMENT}_clu" "${TMP_PATH}/strualn_expanded" ${MERGERESULTBYSET_PAR} \
                || fail "Expand died"
        fi
        INTERMEDIATE="${TMP_PATH}/strualn_expanded"
    fi

    if notExists "${TMP_PATH}/aln.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" tmalign "${QUERY_ALIGNMENT}" "${TARGET_ALIGNMENT}${INDEXEXT}" "${INTERMEDIATE}" "${TMP_PATH}/aln" ${ALIGNMENT_PAR} \
            || fail "Alignment step died"
    fi

    if [ -n "$REMOVE_TMP" ]; then
        echo "Removing temporary files"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/strualn" ${VERBOSITY}
    fi
else

   # 2. Alignment
    if notExists "${TMP_PATH}/strualn.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" $ALIGNMENT_ALGO "${QUERY_ALIGNMENT}" "${TARGET_ALIGNMENT}${INDEXEXT}" "${TMP_PATH}/pref" "${TMP_PATH}/strualn" ${ALIGNMENT_PAR} \
            || fail "Structure alignment step died"
    fi

    if [ -n "${EXPAND}" ]; then
        if notExists "${TMP_PATH}/strualn_expanded.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" mergeresultsbyset "${TMP_PATH}/strualn" "${TARGET_ALIGNMENT}${INDEXEXT}" "${TMP_PATH}/strualn_expanded" ${MERGERESULTBYSET_PAR} \
                || fail "Expand died"
            "$MMSEQS" setextendeddbtype "${TMP_PATH}/strualn_expanded" --extended-dbtype 2
        fi
        INTERMEDIATE="${TMP_PATH}/strualn_expanded"
        if notExists "${TMP_PATH}/aln.dbtype"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" $ALIGNMENT_ALGO "${QUERY_ALIGNMENT}" "${TARGET_ALIGNMENT}${INDEXEXT}" "${TMP_PATH}/strualn_expanded" "${TMP_PATH}/aln" ${ALIGNMENT_PAR} \
                || fail "Alignment step died"
        fi
    else
      "$MMSEQS" mvdb "${TMP_PATH}/strualn" "${TMP_PATH}/aln"
    fi
fi

# shellcheck disable=SC2086
"$MMSEQS" mvdb "${TMP_PATH}/aln" "${RESULTS}" ${VERBOSITY}

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    if [ "$ALIGNMENT_ALGO" = "tmalign" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/strualn" ${VERBOSITY}
    fi

    if [ -n "${EXPAND}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/strualn_expanded" ${VERBOSITY}
    fi

    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/pref" ${VERBOSITY}
fi
