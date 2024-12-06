#!/bin/sh -e

LIST="$1"
IN="$2"
OUT="$3"

if [ ! -d "${TMP_PATH}" ]; then
    mkdir -p "${TMP_PATH}"
fi

if notExists "${OUT}_ca.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "${LIST}" "${IN}" "${OUT}" ${CREATESTRUCTSUBDB_PAR} \
        || fail "createsubdb died"
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "${LIST}" "${IN}_ss" "${OUT}_ss" ${CREATESTRUCTSUBDB_PAR} \
        || fail "createsubdb died"
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "${LIST}" "${IN}_ca" "${OUT}_ca" ${CREATESTRUCTSUBDB_PAR} \
        || fail "createsubdb died"
fi

rm -f "${TMP_PATH}/createstructsubdb.sh"