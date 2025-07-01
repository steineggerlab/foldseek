#!/bin/sh -e

notExists() {
	[ ! -f "$1" ]
}

if [ "${NEEDSET}" -eq "1" ]; then
    # shellcheck disable=SC2086
    $MMSEQS createsimpledb "${INPUT}" "${ALN}_multimerdb" \
        || fail "createsimpledb died"
    INPUT="${ALN}_multimerdb"
fi

if notExists "${RESULT}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" base:clust "${INPUT}" "${ALN}" "${RESULT}" ${CLUST_PAR} \
        || fail "clust died"

    # shellcheck disable=SC2086
    "$MMSEQS" setextendeddbtype "${RESULT}" --extended-dbtype 16 ${VERBOSITY_PAR} \
        || fail "setextendeddbtype died"
fi


