#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

TARGETDB="${INPUT}"
if [ -n "${GPU}" ]; then
    if [ -n "${NOTPADDED}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" makepaddedseqdb "${TARGETDB}" "${TMP_PATH}/target_pad" ${MAKEPADDEDSEQDB_PAR} \
                || fail "makepaddedseqdb died"
        TARGETDB="${TMP_PATH}/target_pad"
    fi
fi

if notExists "${TMP_PATH}/result.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" search "${QUERYDB}" "${TARGETDB}" "${TMP_PATH}/result" "${TMP_PATH}/search_tmp" ${SEARCH_PAR} \
        || fail "Search died"
fi

RESULT="${TMP_PATH}/result"
if [ "$PREFMODE" != "EXHAUSTIVE" ]; then
    if notExists "${TMP_PATH}/result_expand_pref.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" expandmultimer "${QUERYDB}" "${TARGETDB}" "${RESULT}" "${RESULT}_expand_pref" ${THREADS_PAR} \
            || fail "expandmultimer died"
    fi
    if notExists "${TMP_PATH}/result_expand_aligned.dbtype"; then
        if [ "$MULTIMER_ALIGNMENT_ALGO" = "tmalign" ]; then
#            # shellcheck disable=SC2086
#            "$MMSEQS" structurealign "${QUERYDB}" "${TARGETDB}" "${RESULT}_expand_pref" "${RESULT}_expand_aligned_tmp" ${MULTIMER_ALIGN_PREF_PAR} \
#                || fail $MULTIMER_ALIGNMENT_ALGO "died"
            # shellcheck disable=SC2086
            "$MMSEQS" structurealign "${QUERYDB}" "${TARGETDB}" "${RESULT}_expand_pref" "${RESULT}_expand_aligned_tmp" -e 100 ${THREADS_PAR} \
                || fail $MULTIMER_ALIGNMENT_ALGO "died"
            # shellcheck disable=SC2086
            "$MMSEQS" tmalign "${QUERYDB}" "${TARGETDB}" "${RESULT}_expand_aligned_tmp" "${RESULT}_expand_aligned" ${MULTIMER_ALIGN_PAR} \
                || fail $MULTIMER_ALIGNMENT_ALGO "died"
        else
            # shellcheck disable=SC2086
            "$MMSEQS" $MULTIMER_ALIGNMENT_ALGO "${QUERYDB}" "${TARGETDB}" "${RESULT}_expand_pref" "${RESULT}_expand_aligned" ${MULTIMER_ALIGN_PAR} \
                || fail $MULTIMER_ALIGNMENT_ALGO "died"
        fi
    fi
    RESULT="${TMP_PATH}/result_expand_aligned"
fi

if notExists "${TMP_PATH}/scoremultimer.dbtype"; then
    # shellcheck disable=SC2086
    $MMSEQS scoremultimer "${QUERYDB}" "${TARGETDB}" "${RESULT}" "${TMP_PATH}/scoremultimer" ${SCOREMULTIMER_PAR} \
        || fail "scoremultimer died"
    # shellcheck disable=SC2086
    "$MMSEQS" mvdb "${TMP_PATH}/scoremultimer" "${OUTPUT}" \
        || fail "mvdb died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result" ${VERBOSITY}

    if [ "$PREFMODE" != "EXHAUSTIVE" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/result_expand_aligned" ${VERBOSITY}
    fi
    rm -rf "${TMP_PATH}/search_tmp"
    rm -f "${TMP_PATH}/multimersearch.sh"
fi
