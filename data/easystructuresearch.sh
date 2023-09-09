#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

if notExists "${1}.dbtype"; then
    if notExists "${TMP_PATH}/query.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createdb "$@" "${TMP_PATH}/query" ${CREATEDB_PAR} \
            || fail "query createdb died"
    fi
    QUERY="${TMP_PATH}/query"
else
    QUERY="$1"
fi

if notExists "${TARGET}.dbtype"; then
    if notExists "${TMP_PATH}/target"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createdb "${TARGET}" "${TMP_PATH}/target" ${CREATEDB_PAR} \
            || fail "target createdb died"
    fi
    TARGET="${TMP_PATH}/target"
fi


INTERMEDIATE="${TMP_PATH}/result"
if notExists "${INTERMEDIATE}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" search "${QUERY}" "${TARGET}" "${INTERMEDIATE}" "${TMP_PATH}/search_tmp" ${SEARCH_PAR} \
        || fail "Search died"
fi

if [ -n "${GREEDY_BEST_HITS}" ]; then
    if notExists "${TMP_PATH}/result_best.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" summarizeresult "${TMP_PATH}/result" "${TMP_PATH}/result_best" ${SUMMARIZE_PAR} \
            || fail "Search died"
    fi
    INTERMEDIATE="${TMP_PATH}/result_best"
fi

if notExists "${TMP_PATH}/alis.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" convertalis "${QUERY}" "${TARGET}${INDEXEXT}" "${INTERMEDIATE}" "${RESULTS}" ${CONVERT_PAR} \
        || fail "Convert Alignments died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    if [ -n "${GREEDY_BEST_HITS}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/result_best" ${VERBOSITY}
    fi
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result" ${VERBOSITY}
    if [ -z "${LEAVE_INPUT}" ]; then
        if [ -f "${TMP_PATH}/target.dbtype" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_h" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_ca" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_ss" ${VERBOSITY}
        fi
        if [ -f "${TMP_PATH}/query.dbtype" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/query" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/query_h" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/query_ca" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/query_ss" ${VERBOSITY}
        fi
    fi
    rm -rf "${TMP_PATH}/search_tmp"
    rm -f "${TMP_PATH}/easystructuresearch.sh"
fi
