#!/bin/sh -ex
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
    [ ! -f "$1" ]
}

hasCommand () {
    command -v "$1" >/dev/null 2>&1
}

ARR=""
push_back() {
    # shellcheck disable=SC1003
    CURR="$(printf '%s' "$1" | awk '{ gsub(/'\''/, "'\''\\'\'''\''"); print; }')"
    if [ -z "$ARR" ]; then
        ARR=''\'$CURR\'''
    else
        ARR=$ARR' '\'$CURR\'''
    fi
}

STRATEGY=""
if hasCommand aria2c; then STRATEGY="$STRATEGY ARIA"; fi
if hasCommand curl;   then STRATEGY="$STRATEGY CURL"; fi
if hasCommand wget;   then STRATEGY="$STRATEGY WGET"; fi
if [ "$STRATEGY" = "" ]; then
    fail "No download tool found in PATH. Please install aria2c, curl or wget."
fi

downloadFile() {
    URL="$1"
    OUTPUT="$2"
    set +e
    for i in $STRATEGY; do
        case "$i" in
        ARIA)
            FILENAME=$(basename "${OUTPUT}")
            DIR=$(dirname "${OUTPUT}")
            aria2c --max-connection-per-server="$ARIA_NUM_CONN" --allow-overwrite=true -o "$FILENAME" -d "$DIR" "$URL" && return 0
            ;;
        CURL)
            curl -L -o "$OUTPUT" "$URL" && return 0
            ;;
        WGET)
            wget -O "$OUTPUT" "$URL" && return 0
            ;;
        esac
    done
    set -e
    fail "Could not download $URL to $OUTPUT"
}

# check number of input variables
[ "$#" -ne 3 ] && echo "Please provide <selection> <outDB> <tmp>" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && mkdir -p "$3";

SELECTION="$1"
OUTDB="$2"
TMP_PATH="$3"

INPUT_TYPE=""
case "${SELECTION}" in
    "AlphafoldDb")
        if notExists "${TMP_PATH}/alphafolddb.tar.gz"; then
            date "+%s" > "${TMP_PATH}/version"
            downloadFile "http://wwwuser.gwdg.de/~compbiol/foldseek/alphafolddb.tar.gz" "${TMP_PATH}/alphafolddb.tar.gz"
        fi
        tar xvfz "${TMP_PATH}/alphafolddb.tar.gz" -C "${TMP_PATH}"
        push_back "${TMP_PATH}/alphafolddb"
        INPUT_TYPE="FOLDSEEK_DB"
    ;;
    "PDB")
        if notExists "${TMP_PATH}/pdb_seqres.txt.gz"; then
            date "+%s" > "${TMP_PATH}/version"
            downloadFile "https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz" "${TMP_PATH}/pdb_seqres.txt.gz"
        fi
        push_back "${TMP_PATH}/pdb_seqres.txt.gz"
        INPUT_TYPE="FOLDSEEK_DB"
    ;;
esac


if notExists "${OUTDB}.dbtype"; then
case "${INPUT_TYPE}" in
    "FOLDSEEK_DB")
        eval "set -- $ARR"
        IN="${*}"
        # shellcheck disable=SC2086
        "${MMSEQS}" mvdb "${IN}" "${OUTDB}" || fail "mv died"
        # shellcheck disable=SC2086
        "${MMSEQS}" mvdb "${IN}_ss" "${OUTDB}_ss" || fail "mv died"
        # shellcheck disable=SC2086
        "${MMSEQS}" mvdb "${IN}_h" "${OUTDB}_h" || fail "mv died"
        # shellcheck disable=SC2086
        "${MMSEQS}" mvdb "${IN}_ca" "${OUTDB}_ca" || fail "mv died"
    ;;
esac
fi




if [ -n "${TAXONOMY}" ] && notExists "${OUTDB}_mapping"; then
    case "${SELECTION}" in
     *)
       # shellcheck disable=SC2086
       "${MMSEQS}" prefixid "${OUTDB}_h" "${TMP_PATH}/header_pref.tsv" --tsv ${THREADS_PAR} \
           || fail "prefixid died"
       awk '{ match($0, / OX=[0-9]+ /); if (RLENGTH != -1) { print $1"\t"substr($0, RSTART+4, RLENGTH-5); next; } match($0, / TaxID=[0-9]+ /); print $1"\t"substr($0, RSTART+7, RLENGTH-8); }' "${TMP_PATH}/header_pref.tsv" \
           | LC_ALL=C sort -n > "${OUTDB}_mapping"
       rm -f "${TMP_PATH}/header_pref.tsv"
       # shellcheck disable=SC2086
       "${MMSEQS}" createtaxdb "${OUTDB}" "${TMP_PATH}/taxonomy" ${THREADS_PAR} \
           || fail "createtaxdb died"
       ;;
     esac
fi

if notExists "${OUTDB}.version"; then
    mv -f "${TMP_PATH}/version" "${OUTDB}.version"
fi

if [ -n "${REMOVE_TMP}" ]; then
    rm -f "${TMP_PATH}/download.sh"
fi
