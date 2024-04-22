#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

if notExists "${QUERY}.dbtype"; then
    if notExists "${TMP_PATH}/query"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createdb "${QUERY}" "${TMP_PATH}/query" ${CREATEDB_PAR} \
            || fail "query createdb died"
    fi
    QUERY="${TMP_PATH}/query"
fi

if notExists "${TARGET}.dbtype"; then
    if notExists "${TMP_PATH}/target"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createdb "${TARGET}" "${TMP_PATH}/target" ${CREATEDB_PAR} \
            || fail "target createdb died"
    fi
    TARGET="${TMP_PATH}/target"
fi

if notExists "${TMP_PATH}/complex_result.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" multimersearch "${QUERY}" "${TARGET}" "${TMP_PATH}/complex_result" "${TMP_PATH}/complexsearch_tmp" ${COMPLEXSEARCH_PAR} \
    || fail "multimerSearch died"
fi

# shellcheck disable=SC2086
"$MMSEQS" convertalis "${QUERY}" "${TARGET}" "${TMP_PATH}/complex_result" "${OUTPUT}" ${CONVERT_PAR} \
    || fail "Convert Alignments died"

if [ -z "${NO_REPORT}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" createmultimerreport "${QUERY}" "${TARGET}" "${TMP_PATH}/complex_result" "${OUTPUT}_report" ${REPORT_PAR} \
        || fail "createmultimerreport died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result" ${VERBOSITY}
    if [ "$PREFMODE" != "EXHAUSTIVE" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/result_expand_aligned" ${VERBOSITY}
    fi
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/complex_result" ${VERBOSITY}
    if [ -z "${LEAVE_INPUT}" ]; then
        if [ -f "${TMP_PATH}/target" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_h" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_ca" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_ss" ${VERBOSITY}
        fi
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_h" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_ca" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_ss" ${VERBOSITY}
    fi
    rm -rf "${TMP_PATH}/complexsearch_tmp"
    rm -f "${TMP_PATH}/easymultimersearch.sh"
fi