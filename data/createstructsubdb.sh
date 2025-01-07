#!/bin/sh -e

if [ -e "${IN}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "${LIST}" "${IN}" "${OUT}" ${CREATESTRUCTSUBDB_PAR} \
        || fail "createsubdb died"
fi

if [ -e "${IN}_ss.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "${LIST}" "${IN}_ss" "${OUT}_ss" ${CREATESTRUCTSUBDB_PAR} \
        || fail "createsubdb died"
fi

if [ -e "${IN}_ca.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "${LIST}" "${IN}_ca" "${OUT}_ca" ${CREATESTRUCTSUBDB_PAR} \
        || fail "createsubdb died"
fi

# if [ -e "${IN}_h.dbtype" ]; then
#     # shellcheck disable=SC2086
#     "$MMSEQS" base:createsubdb "${LIST}" "${IN}_h" "${OUT}_h" ${CREATESTRUCTSUBDB_PAR} \
#         || fail "createsubdb died"
# fi

if [ -e "${OUT}.sh" ]; then
  rm -f -- "${OUT}.sh"
fi
