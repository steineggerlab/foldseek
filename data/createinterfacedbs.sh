#!/bin/sh -e


if [ -e "${IN}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" createinterfacedb "${IN}" "${OUT}" ${CREATEINTERFACEDB_PAR} \
        || fail "createinterfacedb died"
fi

if [ -e "${IN}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" lndb "${IN}_h" "${OUT}_h"
fi

if [ -e "${IN}.lookup" ]; then
    ln -s "${IN}.lookup" "${OUT}.lookup"
    ln -s "${IN}.source" "${OUT}.source"
fi

if [ -e "${OUT}.sh" ]; then
    rm -f -- "${OUT}.sh"
fi