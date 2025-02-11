#!/bin/sh -e

#getdimeridxlist as idxlist

if [ -e "${IN}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "idxlist" "${IN}" "${OUT}" ${CREATESTRUCTSUBDB_PAR} \
        || fail "createsubdb died"
fi

if [ -e "${IN}_ss.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "idxlist" "${IN}_ss" "${OUT}_ss" ${CREATESTRUCTSUBDB_PAR} \
        || fail "createsubdb died"
fi

if [ -e "${IN}_ca.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "idxlist"  "${IN}_ca" "${OUT}_ca" ${CREATESTRUCTSUBDB_PAR} \
        || fail "createsubdb died"
fi