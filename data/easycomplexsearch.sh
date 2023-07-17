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


SEARCH_RESULT="${TMP_PATH}/result"
if notExists "${SEARCH_RESULT}.dbtype"; then
    # shellcheck disable=SC2086

    "$MMSEQS" search "${QUERY}" "${TARGET}" "${SEARCH_RESULT}" "${TMP_PATH}/search_tmp" ${SEARCH_PAR} \
        || fail "Search died"
fi

SCORECOMPLEX_RESULT="${TMP_PATH}/result2"
if notExists "${SCORECOMPLEX_RESULT}/.dbtype"; then
    # shellcheck disable=SC2086
    $MMSEQS scorecomplex "${QUERY}" "${TARGET}" "${SEARCH_RESULT}" ${SCORECOMPLEX_RESULT} ${SCORECOMPLEX_PAR} \
        || fail "ScoreComplex died"
fi

if notExists "${TMP_PATH}/alis.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" convertalis "${QUERY}" "${TARGET}" "${SCORECOMPLEX_RESULT}" "${OUTPUT}" ${CONVERT_PAR} \
        || fail "Convert Alignments died"
fi

## shellcheck disable=SC2086
#"$MMSEQS" createcomplexreport "${QUERY}" "${TARGET}" "${SCORECOMPLEX_RESULT}" "${OUTPUT}" ${CONVERT_PAR} \
#    || fail "Convert Alignments died"






if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result" ${VERBOSITY}
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
    rm -rf "${TMP_PATH}/search_tmp"
    rm -f "${TMP_PATH}/easyscorecomplex.sh"
fi
