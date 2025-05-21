#!/bin/sh -e

if [ -e "${IN}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" createstructsimpledb "${IN}" "${OUT}" ${CREATESIMPLEDB_PAR} --dbtype 0 \
        || fail "createstructsimpledb died"
fi

if [ -e "${IN}_ss.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" createstructsimpledb "${IN}_ss" "${OUT}_ss" ${CREATESIMPLEDB_PAR} --dbtype 0 \
        || fail "createstructsimpledb died"
fi

if [ -e "${IN}_ca.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" createstructsimpledb "${IN}_ca" "${OUT}_ca" ${CREATESIMPLEDB_PAR} --dbtype 101 \
        || fail "createstructsimpledb died"
fi

# make _h, .lookup, .source