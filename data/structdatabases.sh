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

downloadFileList() {
    URL="$1"
    OUTPUT_DIR="$2"
    INPUT_FILE="$OUTPUT_DIR/input.txt"
    downloadFile "$URL" "$INPUT_FILE"
    set +e
    for i in $STRATEGY; do
        case "$i" in
        ARIA)
            aria2c -c --max-connection-per-server="$ARIA_NUM_CONN" --allow-overwrite=true --dir="$OUTPUT_DIR" --input-file="$INPUT_FILE" && return 0
            ;;
        CURL)
            (cd "$OUTPUT_DIR"; xargs -n 1 curl -C - -L -O < "$INPUT_FILE") && return 0
            ;;
        WGET)
            wget --continue -P "$OUTPUT_DIR" --input-file="$INPUT_FILE" && return 0
            ;;
        esac
    done
    set -e
    rm -f "$OUTPUT/input.txt"
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
    "Alphafold/UniProt")
        if notExists "${TMP_PATH}/afdb.tar.gz"; then
            downloadFile "https://foldseek.steineggerlab.workers.dev/afdb.tar.gz" "${TMP_PATH}/afdb.tar.gz"
            downloadFile "https://foldseek.steineggerlab.workers.dev/afdb.version" "${TMP_PATH}/version"
        fi
        tar xvfz "${TMP_PATH}/afdb.tar.gz" -C "${TMP_PATH}"
        push_back "${TMP_PATH}/afdb"
        INPUT_TYPE="FOLDSEEK_DB"
    ;;
    "Alphafold/UniProt50-minimal")
        if notExists "${TMP_PATH}/afdb50.tar.gz"; then
            downloadFile "https://foldseek.steineggerlab.workers.dev/afdb50.tar.gz" "${TMP_PATH}/afdb50.tar.gz"
            downloadFile "https://foldseek.steineggerlab.workers.dev/afdb50.version" "${TMP_PATH}/version"
        fi
        tar xvfz "${TMP_PATH}/afdb50.tar.gz" -C "${TMP_PATH}"
        push_back "${TMP_PATH}/afdb50"
        INPUT_TYPE="FOLDSEEK_DB"
    ;;
    "Alphafold/UniProt50")
        if notExists "${TMP_PATH}/afdb50.tar.gz"; then
            downloadFile "https://foldseek.steineggerlab.workers.dev/afdb50.tar.gz" "${TMP_PATH}/afdb50.tar.gz"
            downloadFile "https://foldseek.steineggerlab.workers.dev/afdb50clusearch.tar.gz" "${TMP_PATH}/afdb50clusearch.tar.gz"
            downloadFile "https://foldseek.steineggerlab.workers.dev/afdb50.version" "${TMP_PATH}/version"
        fi
        tar xvfz "${TMP_PATH}/afdb50.tar.gz" -C "${TMP_PATH}"
        tar xvfz "${TMP_PATH}/afdb50clusearch.tar.gz" -C "${TMP_PATH}"
        push_back "${TMP_PATH}/afdb50"
        INPUT_TYPE="FOLDSEEK_DB"
    ;;
    "Alphafold/Proteome")
        if notExists "${TMP_PATH}/alphafolddb.tar.gz"; then
            downloadFile "https://foldseek.steineggerlab.workers.dev/afdb_proteome.tar.gz" "${TMP_PATH}/afdb_proteome.tar.gz"
            downloadFile "https://foldseek.steineggerlab.workers.dev/afdb_proteome.version" "${TMP_PATH}/version"
        fi
        tar xvfz "${TMP_PATH}/afdb_proteome.tar.gz" -C "${TMP_PATH}"
        push_back "${TMP_PATH}/afdb_proteome"
        INPUT_TYPE="FOLDSEEK_DB"
    ;;
    "Alphafold/Swiss-Prot")
        if notExists "${TMP_PATH}/alphafold_swissprot.tar.gz"; then
            downloadFile "https://foldseek.steineggerlab.workers.dev/afdb_swissprot.tar.gz" "${TMP_PATH}/afdb_swissprot.tar.gz"
            downloadFile "https://foldseek.steineggerlab.workers.dev/afdb_swissprot.version" "${TMP_PATH}/version"
        fi
        tar xvfz "${TMP_PATH}/afdb_swissprot.tar.gz" -C "${TMP_PATH}"
        push_back "${TMP_PATH}/afdb_swissprot"
        INPUT_TYPE="FOLDSEEK_DB"
    ;;
    "ESMAtlas30")
        downloadFileList "https://raw.githubusercontent.com/facebookresearch/esm/main/scripts/atlas/v0/highquality_clust30/foldseekdb.txt" "${TMP_PATH}/"
        printf "v0 %s\n" "$(date "+%s")" > "${TMP_PATH}/version"
        push_back "${TMP_PATH}/highquality_clust30"
        INPUT_TYPE="FOLDSEEK_DB"
    ;;
    "PDB")
        if notExists "${TMP_PATH}/pdb.tar.gz"; then
            downloadFile "https://foldseek.steineggerlab.workers.dev/pdb100.tar.gz" "${TMP_PATH}/pdb.tar.gz"
            downloadFile "https://foldseek.steineggerlab.workers.dev/pdb100.version" "${TMP_PATH}/version"
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
        for SUFFIX in ".source" "_mapping" "_taxonomy"; do
            if [ -e "${IN}_seq${SUFFIX}" ]; then
                mv -f -- "${IN}_seq${SUFFIX}" "${OUTDB}_seq${SUFFIX}"
            fi
            if [ -e "${IN}${SUFFIX}" ]; then
                mv -f -- "${IN}${SUFFIX}" "${OUTDB}${SUFFIX}"
            fi
        done

        for SUFFIX in "" "_ss" "_h" "_ca"; do
            if [ -e "${IN}_seq${SUFFIX}.dbtype" ]; then
                # shellcheck disable=SC2086
                "${MMSEQS}" mvdb "${IN}_seq${SUFFIX}" "${OUTDB}_seq${SUFFIX}" || fail "mv died"
            fi
            # shellcheck disable=SC2086
            "${MMSEQS}" mvdb "${IN}${SUFFIX}" "${OUTDB}${SUFFIX}" || fail "mv died"
        done

        if [ -e "${IN}_clu.dbtype" ]; then
            # shellcheck disable=SC2086
            "${MMSEQS}" mvdb "${IN}_clu" "${OUTDB}_clu" || fail "mv died"
        fi
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
        if [ -e "${IN}_seq.dbtype" ]; then
            mv -f -- "${IN}_seq_mapping" "${OUTDB}_seq_mapping"
            mv -f -- "${IN}_seq_taxonomy" "${OUTDB}_seq_taxonomy"
        fi
    ;;
esac
fi

if notExists "${OUTDB}.version"; then
    mv -f "${TMP_PATH}/version" "${OUTDB}.version"
fi

if [ -n "${REMOVE_TMP}" ]; then
    rm -f "${TMP_PATH}/download.sh"
fi
