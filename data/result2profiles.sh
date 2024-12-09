#!/bin/sh -e

IN1="$1"
IN2="$2"
RESULT="$3"
OUT="$4"

if [ ! -d "${TMP_PATH}" ]; then
    mkdir -p "${TMP_PATH}"
fi

# shellcheck disable=SC2086
"$MMSEQS" base:result2profile "${IN1}" "${IN2}" "${RESULT}" "${OUT}" ${RESULT2PROFILES_PAR} \
    || fail "result2profile died"

# shellcheck disable=SC2086
"$MMSEQS" base:result2profile "${IN1}_ss" "${IN2}_ss" "${RESULT}" "${OUT}_ss" ${RESULT2PROFILES_PAR} \
    || fail "result2profile died"
# shellcheck disable=SC2086
"$MMSEQS" lndb "${IN1}_ca"  "${OUT}_ca" ${VERBOSITY} \
    || fail "Create lndb died"

rm -f "${TMP_PATH}/result2profiles.sh"