#!/bin/sh -e

INPUT="$1"
RESULT="$3"
OUT="$4"

if [ -e "${INPUT}.dbtype" ]; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" base:result2profile "${INPUT}" "${TARGET}${INDEXEXT}" "${RESULT}" "${OUT}" ${PROFILE_PAR} \
        || fail "result2profile died"
fi

if [ -e "${INPUT}_ss.dbtype" ]; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" base:result2profile "${INPUT}_ss" "${TARGET}_ss${INDEXEXT}" "${RESULT}" "${OUT}_ss" ${PROFILE_SS_PAR} \
        || fail "result2profile died"
fi

if [ -e "${INPUT}_ca.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" lndb "${INPUT}_ca" "${OUT}_ca" ${VERBOSITY} \
        || fail "Create lndb died"
fi

if [ -e "${INPUT}_h.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" lndb "${INPUT}_h" "${OUT}_h" ${VERBOSITY} \
        || fail "Create lndb died"
fi

if [ -e "${INPUT}.sh" ]; then
    rm -f -- "${OUT}.sh"
fi
