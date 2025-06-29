#!/bin/sh -e


if [ "${NEEDSET}" -eq "1" ]; then
    # shellcheck disable=SC2086
    $MMSEQS createsimpledb "${INPUT}" "${ALN}_multimerdb" ${VERBOSITY} \
        || fail "createsimpledb died"
    INPUT="${ALN}_multimerdb"
fi

if notExists "${RESULT}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" clust "${INPUT}" "${ALN}" "${RESULT}" ${CLUST_PAR} \
        || fail "clust died"
fi
