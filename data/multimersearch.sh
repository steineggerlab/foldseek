#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

if notExists "${TMP_PATH}/result.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" search "${QUERYDB}" "${TARGETDB}" "${TMP_PATH}/result" "${TMP_PATH}/search_tmp" ${SEARCH_PAR} \
        || fail "Search died"
fi

RESULT="${TMP_PATH}/result"
if [ "$PREFMODE" != "EXHAUSTIVE" ]; then
    if notExists "${TMP_PATH}/result_expand_pref.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" expandcomplex "${QUERYDB}" "${TARGETDB}" "${RESULT}" "${TMP_PATH}/result_expand_pref" ${THREADS_PAR} \
            || fail "Expandcomplex died"
    fi
    if notExists "${TMP_PATH}/result_expand_aligned.dbtype"; then
        if [ "$COMPLEX_ALIGNMENT_ALGO" = "tmalign" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" structurealign "${QUERYDB}" "${TARGETDB}" "${TMP_PATH}/result_expand_pref" "${TMP_PATH}/result_expand_aligned_tmp" -e 100 \
                || fail $COMPLEX_ALIGNMENT_ALGO "died"
            # shellcheck disable=SC2086
            "$MMSEQS" tmalign "${QUERYDB}" "${TARGETDB}" "${TMP_PATH}/result_expand_aligned_tmp" "${TMP_PATH}/result_expand_aligned" ${COMPLEX_ALIGN_PAR} \
                || fail $COMPLEX_ALIGNMENT_ALGO "died"
        else
            # shellcheck disable=SC2086
            "$MMSEQS" $COMPLEX_ALIGNMENT_ALGO "${QUERYDB}" "${TARGETDB}" "${TMP_PATH}/result_expand_pref" "${TMP_PATH}/result_expand_aligned" ${COMPLEX_ALIGN_PAR} \
                || fail $COMPLEX_ALIGNMENT_ALGO "died"
        fi
    fi
    RESULT="${TMP_PATH}/result_expand_aligned"
fi
if notExists "${TMP_PATH}/complex_result.dbtype"; then
    # shellcheck disable=SC2086
    $MMSEQS scorecomplex "${QUERYDB}" "${TARGETDB}" "${RESULT}" "${OUTPUT}" ${SCORECOMPLEX_PAR} \
        || fail "ScoreComplex died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result" ${VERBOSITY}
    if [ "$PREFMODE" != "EXHAUSTIVE" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/result_expand_aligned" ${VERBOSITY}
    fi
    rm -rf "${TMP_PATH}/search_tmp"
    rm -f "${TMP_PATH}/complexsearch.sh"
fi
