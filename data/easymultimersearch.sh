#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

exists() {
	[ -f "$1" ]
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

    if [ -n "${GPU}" ]; then
        if notExists "${TMP_PATH}/target_pad"; then
            # shellcheck disable=SC2086
            "$MMSEQS" makepaddedseqdb "${TMP_PATH}/target" "${TMP_PATH}/target_pad" ${MAKEPADDEDSEQDB_PAR} \
                || fail "makepaddedseqdb died"
        fi
        TARGET="${TMP_PATH}/target_pad"
    fi
fi

if notExists "${TMP_PATH}/multimer_result.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" multimersearch "${QUERY}" "${TARGET}" "${TMP_PATH}/multimer_result" "${TMP_PATH}/multimersearch_tmp" ${MULTIMERSEARCH_PAR} \
    || fail "multimersearch died"
fi

# # shellcheck disable=SC2086
# "$MMSEQS" createtsv "${TMP_PATH}/multimer_result_query_multimerdb" "${TMP_PATH}/multimer_result_target_multimerdb" "${TMP_PATH}/multimer_result" "${OUTPUT}" ${THREADS_PAR} \
#     || fail "createtsv died"
# shellcheck disable=SC2086
"$MMSEQS" createtsv "${QUERY}" "${TARGET}" "${TMP_PATH}/multimer_result" "${OUTPUT}" ${THREADS_PAR} \
    || fail "createtsv died"

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/multimer_result" ${VERBOSITY}
    if [ -z "${LEAVE_INPUT}" ]; then
        if [ -f "${TMP_PATH}/target" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_h" ${VERBOSITY}
            if exists "${TMP_PATH}/target_ca.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/target_ca" ${VERBOSITY}
            fi
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_ss" ${VERBOSITY}
        fi
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_h" ${VERBOSITY}
        if exists "${TMP_PATH}/target_ca.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/query_ca" ${VERBOSITY}
        fi
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_ss" ${VERBOSITY}
    fi
    rm -rf "${TMP_PATH}/multimersearch_tmp"
    rm -f "${TMP_PATH}/easymultimersearch.sh"
fi