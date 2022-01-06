#!/bin/sh -e
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
    "Alphafold/Proteome")
        if notExists "${TMP_PATH}/alphafolddb.tar.gz"; then
            date "+%s" > "${TMP_PATH}/version"
            downloadFile "https://wwwuser.gwdg.de/~compbiol/foldseek/alphafolddb.tar.gz" "${TMP_PATH}/alphafolddb.tar.gz"
        fi
        tar xvfz "${TMP_PATH}/alphafolddb.tar.gz" -C "${TMP_PATH}"
        push_back "${TMP_PATH}/alphafolddb"
        INPUT_TYPE="FOLDSEEK_DB"
    ;;
    "Alphafold/Swiss-Prot")
          if notExists "${TMP_PATH}/alphafold_swissprot.tar.gz"; then
              date "+%s" > "${TMP_PATH}/version"
              downloadFile "https://wwwuser.gwdg.de/~compbiol/foldseek/alphafold_swissprot.tar.gz" "${TMP_PATH}/alphafold_swissprot.tar.gz"
          fi
          tar xvfz "${TMP_PATH}/alphafold_swissprot.tar.gz" -C "${TMP_PATH}"
          push_back "${TMP_PATH}/alphafold_swissprot"
          INPUT_TYPE="FOLDSEEK_DB"
    ;;
    "PDB")
        if notExists "${TMP_PATH}/pdb.tar.gz"; then
            date "+%s" > "${TMP_PATH}/version"
            downloadFile "https://wwwuser.gwdg.de/~compbiol/foldseek/pdb.tar.gz" "${TMP_PATH}/pdb.tar.gz"
        fi
        tar xvfz "${TMP_PATH}/pdb.tar.gz" -C "${TMP_PATH}"
        push_back "${TMP_PATH}/pdb"
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
case "${INPUT_TYPE}" in
    "FOLDSEEK_DB")
        eval "set -- $ARR"
        IN="${*}"
        mv -f -- "${IN}_mapping" "${OUTDB}_mapping"
        mv -f -- "${IN}_taxonomy" "${OUTDB}_taxonomy"
    ;;
esac
fi

if notExists "${OUTDB}.version"; then
    mv -f "${TMP_PATH}/version" "${OUTDB}.version"
fi

if [ -n "${REMOVE_TMP}" ]; then
    rm -f "${TMP_PATH}/download.sh"
fi
