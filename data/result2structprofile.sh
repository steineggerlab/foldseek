#!/bin/sh -e

IN1="$1"
IN2="$2"
RESULT="$3"
OUT="$4"

if [ -e "${IN1}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" base:result2profile "${IN1}" "${IN2}" "${RESULT}" "${OUT}" ${PROFILE_PAR} \
        || fail "result2profile died"
fi

if [ -e "${IN1}_ss.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" base:result2profile "${IN1}_ss" "${IN2}_ss" "${RESULT}" "${OUT}_ss" ${PROFILE_SS_PAR} \
        || fail "result2profile died"
fi

if [ -e "${IN1}_ca.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" lndb "${IN1}_ca"  "${OUT}_ca" ${VERBOSITY} \
        || fail "Create lndb died"
fi

if [ -e "${IN1}_h.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" lndb "${IN1}_h"  "${OUT}_h" ${VERBOSITY} \
        || fail "Create lndb died"
fi

if [ -e "${OUT}.sh" ]; then
  rm -f -- "${OUT}.sh"
fi
