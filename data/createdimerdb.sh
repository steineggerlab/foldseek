#!/bin/sh -e


#currently checking if shraing monomers among dimers is possible. (So the lookup file's first column might not be unique)
#it seems working. I'll update this code later


# TODO: getmonmer index list as idxlist

if [ -e "${IN}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "idxlist" "${IN}" "${OUT}" --subdb-mode 0 ${VERBOSITY_PAR} \
        || fail "createsubdb died"
    rm -- "${OUT}_h" "${OUT}_h.dbtype" "${OUT}_h.index"
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "idxlist" "${IN}_h" "${OUT}_h" --subdb-mode 0  ${VERBOSITY_PAR} \
        || fail "createsubdb died"
    rm -- "${OUT}.lookup" "${OUT}.source" 
fi

# TODO: make source file of all dimers
# TODO: make lookup file reflecting all monomer-dimer relation
