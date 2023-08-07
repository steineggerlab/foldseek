#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

DB="$1"
TMP_PATH="$2"

# shellcheck disable=SC2086
"$MMSEQS" mmcreateindex "${DB}" "${TMP_PATH}" ${CREATEINDEX_PAR} --index-subset 2 \
    || fail "createindex died"

if notExists "${DB}_ss_h.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" lndb "${DB}_h" "${DB}_ss_h" ${VERBOSITY_PAR} \
        || fail "lndb header died"
fi

IS_CLUDDB="0"
if [ -f "${DB}_clu.dbtype" ]; then
   IS_CLUDDB="1"
    # shellcheck disable=SC2086
    "$MMSEQS" lndb "${DB}_clu" "${DB}_ss_clu" ${VERBOSITY_PAR} \
        || fail "lndb ss_clu died"
fi

if [ -f "${DB}_aln.dbtype" ]; then
    IS_CLUDDB="1"
    # shellcheck disable=SC2086
    "$MMSEQS" lndb "${DB}_aln" "${DB}_ss_aln" ${VERBOSITY_PAR} \
        || fail "lndb ss_aln died"
fi

# shellcheck disable=SC2086
"$MMSEQS" mmcreateindex "${DB}_ss" "${TMP_PATH}" ${CREATEINDEX_PAR} --index-subset ${SS_SUBSET_MODE}  --index-dbsuffix "_ss" \
    || fail "createindex died"

if [ -n "$INCLUDE_CA" ]; then
  # shellcheck disable=SC2086
  "$MMSEQS" appenddbtoindex "${DB}_ca" "${DB}.idx" --id-list ${INDEX_DB_CA_KEY_DB1} ${VERBOSITY_PAR} \
      || fail "appenddbtoindex died"
  if [ -f "${DB}_seq_ca.dbtype" ] && [ "$IS_CLUDDB" = "1" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" appenddbtoindex "${DB}_seq_ca" "${DB}.idx" --id-list ${INDEX_DB_CA_KEY_DB2} ${VERBOSITY_PAR} \
        || fail "appenddbtoindex died"
  fi
fi