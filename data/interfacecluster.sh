#!/bin/sh -e

exists() {
	[ -f "$1" ]
}

fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

QUERYDB="${INPUT}"
if exists "${INPUT}.dbtype" ; then
    if [ -n "${ISINTERFACEDB}" ]; then
        if [ -n "${GPU}" ]; then
            if [ -n "${NOTPADDED}" ]; then 
                # shellcheck disable=SC2086
                "$MMSEQS" makepaddedseqdb "${QUERYDB}" "${TMP_PATH}/input_pad" ${MAKEPADDEDSEQDB_PAR} \
                        || fail "makepaddedseqdb died"
                QUERYDB="${TMP_PATH}/input_pad" 
            fi
        fi
    else
        if [ -z "${NOTPADDED}" ]; then
            fail "We cannot make an interface db out of padded db"
        fi
        # shellcheck disable=SC2086
        "$MMSEQS" createdimerdb "${QUERYDB}" "${TMP_PATH}/dimerdb" "${TMP_PATH}/dimertmp" ${THREADS_PAR} \
            || fail "createdimerdb died"
        # shellcheck disable=SC2086
        "$MMSEQS" createinterfacedb "${TMP_PATH}/dimerdb" "${TMP_PATH}/interfacedb" ${THREADS_PAR} \
            || fail "createinterfacedb died"
        QUERYDB="${TMP_PATH}/interfacedb"
        if [ -n "${GPU}" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" makepaddedseqdb "${TMP_PATH}/interfacedb" "${TMP_PATH}/interfacedb_pad" ${MAKEPADDEDSEQDB_PAR} \
                    || fail "makepaddedseqdb died"
            QUERYDB="${TMP_PATH}/interfacedb_pad"
        fi
    fi
    # shellcheck disable=SC2086
    "$MMSEQS" multimercluster "${QUERYDB}" "${OUT}" "${TMP_PATH}/multimercluster_tmp" ${INTERFACECLUSTER_PAR} \
        || fail "multimercluster died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    if [ -z "${ISINTERFACEDB}" ]; then
        if [ -n "${GPU}" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/interfacedb_pad" ${VERBOSITY_PAR} 
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/interfacedb_pad_ss" ${VERBOSITY_PAR} 
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/interfacedb_pad_ca" ${VERBOSITY_PAR} 
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/interfacedb_pad_h" ${VERBOSITY_PAR} 
        fi
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/interfacedb" ${VERBOSITY_PAR} 
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/interfacedb_h" ${VERBOSITY_PAR} 
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/interfacedb_ss" ${VERBOSITY_PAR} 
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/interfacedb_ca" ${VERBOSITY_PAR} 
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/dimerdb" ${VERBOSITY_PAR} 
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/dimerdb_h" ${VERBOSITY_PAR} 
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/dimerdb_ss" ${VERBOSITY_PAR} 
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/dimerdb_ca" ${VERBOSITY_PAR}
        rm -rf "${TMP_PATH}/dimertmp"
    fi
    rm -rf -- "${TMP_PATH}/multimercluster_tmp"
    rm -f -- "${TMP_PATH}/interfacecluster.sh"
fi