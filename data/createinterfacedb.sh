#!/bin/sh -e

exists() {
	[ -f "$1" ]
}

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [ -z "${1##*/*}" ]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

if [ -e "${IN}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" createStructinterfacedb "${IN}" "${OUT}" ${CREATESTRUCTINTERFACEDB_PAR} \
        || fail "createStructinterfacedb died"
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${OUT}.index" "${IN}_h" "${OUT}_h" \
        || fail "createsubdb died"
    ln -sf "$(abspath "${IN}.lookup")" "${OUT}.lookup"
    ln -sf "$(abspath "${IN}.source")" "${OUT}.source"
    if exists "${IN}_mapping"; then
        ln -sf "$(abspath "${IN}_mapping")" "${OUT}_mapping"
    fi
    if exists "${IN}_taxonomy"; then
        ln -sf "$(abspath "${IN}_taxonomy")" "${OUT}_taxonomy"
    fi
fi

if [ -e "${OUT}.sh" ]; then
    rm -f -- "${OUT}.sh"
fi