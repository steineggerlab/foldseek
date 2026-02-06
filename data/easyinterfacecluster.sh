#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

exists() {
	[ -f "$1" ]
}

INTERFACEDB="${INPUT}"
if exists "${INPUT}.dbtype"; then
    if [ -n "${ISINTERFACEDB}" ]; then
        if [ -n "${GPU}" ]; then
            if [ -n "${NOTPADDED}" ]; then
                # shellcheck disable=SC2086
                "$MMSEQS" makepaddedseqdb "${INTERFACEDB}" "${TMP_PATH}/interfacedb_pad" ${MAKEPADDEDSEQDB_PAR} \
                    || fail "makepaddedseqdb died"
                INTERFACEDB="${TMP_PATH}/interfacedb_pad"
            fi
        fi
    else
        # shellcheck disable=SC2086
        "$MMSEQS" createdimerdb "${INTERFACEDB}" "${TMP_PATH}/dimerdb" "${TMP_PATH}/dimertmp" ${THREADS_PAR} \
            || fail "createdimerdb died"
        # shellcheck disable=SC2086
        "$MMSEQS" createinterfacedb "${TMP_PATH}/dimerdb" "${TMP_PATH}/interfacedb" ${THREADS_PAR} \
            || fail "createinterfacedb died"
        INTERFACEDB="${TMP_PATH}/interfacedb"
        if [ -n "${GPU}" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" makepaddedseqdb "${INTERFACEDB}" "${TMP_PATH}/interfacedb_pad" ${MAKEPADDEDSEQDB_PAR} \
                || fail "makepaddedseqdb died"
            INTERFACEDB="${TMP_PATH}/interfacedb_pad"
        fi
    fi
else
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "${INPUT}" "${TMP_PATH}/query" ${CREATEDB_PAR} \
        || fail "query createdb died"
    # shellcheck disable=SC2086
    "$MMSEQS" createdimerdb "${TMP_PATH}/query" "${TMP_PATH}/dimerdb" "${TMP_PATH}/dimertmp" ${THREADS_PAR} \
        || fail "createdimerdb died"
    # shellcheck disable=SC2086
    "$MMSEQS" createinterfacedb "${TMP_PATH}/dimerdb" "${TMP_PATH}/interfacedb" ${THREADS_PAR} \
        || fail "createinterfacedb died"
    INTERFACEDB="${TMP_PATH}/interfacedb"
    if [ -n "${GPU}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" makepaddedseqdb "${INTERFACEDB}" "${TMP_PATH}/interfacedb_pad" ${MAKEPADDEDSEQDB_PAR} \
            || fail "makepaddedseqdb died"
        INTERFACEDB="${TMP_PATH}/interfacedb_pad"
    fi
fi

# shellcheck disable=SC2086
"$MMSEQS" easy-multimercluster "${INTERFACEDB}" "${TMP_PATH}/interface_clu" "${TMP_PATH}/interfacecluster_tmp" ${EASYMULTIMERCLUSTER_PAR} \
    || fail "Interfacecluster died"
mv -f -- "${TMP_PATH}/interface_clu_rep_seq.fasta" "${RESULT}_rep_seq.fasta"
mv -f -- "${TMP_PATH}/interface_clu_cluster.tsv" "${RESULT}_cluster.tsv"

if [ -n "${REMOVE_TMP}" ]; then
    if notExists "${INPUT}.dbtype"; then
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
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query" ${VERBOSITY_PAR} 
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_h" ${VERBOSITY_PAR} 
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_ca" ${VERBOSITY_PAR} 
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_ss" ${VERBOSITY_PAR} 
        rm -rf -- "${TMP_PATH}/dimertmp"
    elif [ -z "${ISINTERFACEDB}" ]; then
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
        rm -rf -- "${TMP_PATH}/dimertmp"
    fi
    rm -rf -- "${TMP_PATH}/interfacecluster_tmp"
    rm -r -- "${TMP_PATH}/easyinterfacecluster.sh"
fi
