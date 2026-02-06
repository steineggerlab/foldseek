#!/bin/sh -e

exists() {
	[ -f "$1" ]
}

fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

TARGET="${INPUT_TARGET}"
QUERY="${INPUT_QUERY}"
if exists "${TARGET}.dbtype" ; then
    if [ -n "${ISINTERFACEDB_TARGET}" ]; then
        if [ -n "${GPU}" ]; then
            if [ -n "${NOTPADDED}" ]; then 
                # shellcheck disable=SC2086
                "$MMSEQS" makepaddedseqdb "${TARGET}" "${TMP_PATH}/interfacedb_target_pad" ${MAKEPADDEDSEQDB_PAR} \
                        || fail "makepaddedseqdb died"
                TARGET="${TMP_PATH}/interfacedb_target_pad" 
            fi
        fi
    else
        if [ -z "${NOTPADDED_TARGET}" ]; then
            fail "We cannot make an interface db out of padded db"
        fi
        # shellcheck disable=SC2086
        "$MMSEQS" createdimerdb "${TARGET}" "${TMP_PATH}/dimerdb_target" "${TMP_PATH}/dimertmp_target" ${THREADS_PAR} \
            || fail "createdimerdb died"
        # shellcheck disable=SC2086
        "$MMSEQS" createinterfacedb "${TMP_PATH}/dimerdb_target" "${TMP_PATH}/interfacedb_target" ${THREADS_PAR} \
            || fail "createinterfacedb died"
        TARGET="${TMP_PATH}/interfacedb_target"
        if [ -n "${GPU}" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" makepaddedseqdb "${TMP_PATH}/interfacedb_target" "${TMP_PATH}/interfacedb_target_pad" ${MAKEPADDEDSEQDB_PAR} \
                    || fail "makepaddedseqdb died"
            TARGET="${TMP_PATH}/interfacedb_target_pad"
        fi
    fi
else
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "${TARGET}" "${TMP_PATH}/target" ${CREATEDB_PAR} \
        || fail "target createdb died"
    # shellcheck disable=SC2086
    "$MMSEQS" createdimerdb "${TMP_PATH}/target" "${TMP_PATH}/dimerdb_target" "${TMP_PATH}/dimertmp_target" ${THREADS_PAR} \
        || fail "createdimerdb died"
    # shellcheck disable=SC2086
    "$MMSEQS" createinterfacedb "${TMP_PATH}/dimerdb_target" "${TMP_PATH}/interfacedb_target" ${THREADS_PAR} \
        || fail "createinterfacedb died"
    TARGET="${TMP_PATH}/interfacedb_target"
    if [ -n "${GPU}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" makepaddedseqdb "${TARGET}" "${TMP_PATH}/interfacedb_target_pad" ${MAKEPADDEDSEQDB_PAR} \
            || fail "makepaddedseqdb died"
        TARGET="${TMP_PATH}/interfacedb_target_pad"
    fi
fi
if exists "${QUERY}.dbtype" ; then
    if [ -z "${ISINTERFACEDB_QUERY}" ]; then
        if [ -z "${NOTPADDED_QUERY}" ]; then
            fail "We cannot make an interface db out of padded db"
        fi
        # shellcheck disable=SC2086
        "$MMSEQS" createdimerdb "${QUERY}" "${TMP_PATH}/dimerdb_query" "${TMP_PATH}/dimertmp_query" ${THREADS_PAR} \
            || fail "createdimerdb died"
        # shellcheck disable=SC2086
        "$MMSEQS" createinterfacedb "${TMP_PATH}/dimerdb_query" "${TMP_PATH}/interfacedb_query" ${THREADS_PAR} \
            || fail "createinterfacedb died"
        QUERY="${TMP_PATH}/interfacedb_query"
    fi
else
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "${QUERY}" "${TMP_PATH}/query" ${CREATEDB_PAR} \
        || fail "query createdb died"
    # shellcheck disable=SC2086
    "$MMSEQS" createdimerdb "${TMP_PATH}/query" "${TMP_PATH}/dimerdb_query" "${TMP_PATH}/dimertmp_query" ${THREADS_PAR} \
        || fail "createdimerdb died"
    # shellcheck disable=SC2086
    "$MMSEQS" createinterfacedb "${TMP_PATH}/dimerdb_query" "${TMP_PATH}/interfacedb_query" ${THREADS_PAR} \
        || fail "createinterfacedb died"
    QUERY="${TMP_PATH}/interfacedb_query"
fi

# shellcheck disable=SC2086
"$MMSEQS" easy-multimersearch "${QUERY}" "${TARGET}" "${OUT}" "${TMP_PATH}/easymultimersearch_tmp" ${EASYMULTIMERSEARCH_PAR} \
    || fail "multimersearch died"

if [ -n "${REMOVE_TMP}" ]; then
    # TODO: remove all the dbs
    rm -rf -- "${TMP_PATH}/easymultimersearch_tmp"
    rm -f -- "${TMP_PATH}/easyinterfacesearch.sh"
fi