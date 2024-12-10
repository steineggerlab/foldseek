#!/bin/sh -e

IN="$1"
OUT="$2"

outname=$(basename "${OUT}")
gpu_mapping="${outname}_ss.gpu_mapping"

if [ ! -d "${TMP_PATH}" ]; then
    mkdir -p "${TMP_PATH}"
fi

# shellcheck disable=SC2086
"$MMSEQS" lndb "${IN}_h"  "${IN}_ss_h" ${VERBOSITY} \
    || fail "lndb died"

# shellcheck disable=SC2086
"$MMSEQS" base:makepaddedseqdb "${IN}_ss" "${OUT}_ss" ${MAKEPADDEDSEQDB_PAR} \
    || fail "mmseqs makepaddedseqdb died"

awk '{ print $3"\t"$1; }' "${OUT}_ss.lookup" > "${TMP_PATH}/${gpu_mapping}"

# shellcheck disable=SC2086
"$MMSEQS" renamedbkeys "${TMP_PATH}/${gpu_mapping}" "${IN}" "${OUT}" \
    --subdb-mode 1 ${THREADS_PAR} \
    || fail "renamedbkeys died"

# shellcheck disable=SC2086
"$MMSEQS" renamedbkeys "${TMP_PATH}/${gpu_mapping}" "${IN}_ca" "${OUT}_ca" \
    --subdb-mode 1 ${THREADS_PAR} \
    || fail "renamedbkeys died"

rm -f -- "${TMP_PATH}/${gpu_mapping}"
rm -f "${TMP_PATH}/makepaddeddb.sh"