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

mapCmplName2ChainKeys() {
    awk -F"\t" 'FNR==1 {++fIndex}
        fIndex==1 {
            repName[$1]=1
            if (match($1, /MODEL/)){
                tmpName[$1]=1
            }else{
                tmpName[$1"_MODEL_1"]=1 
            }
            next
        }
        fIndex==2{
            if (match($2, /MODEL/)){
                if ($2 in tmpName){
                repId[$1]=1
                }else{
                    ho[1]=1
                }
            }else{
                if ($2 in repName){
                repId[$1]=1
                }
            }
            next
        }
        {
            if ($3 in repId){
                print $1
            }
        }
    ' "${1}" "${2}.source" "${2}.lookup" > "${3}"
}

checkifIndb() {
    awk 'FNR == NR { name[$1] = 1; next } $1 in name { print }' "${1}.index" "${2}" > "${3}"
}

postprocessFasta() {
    awk ' BEGIN {FS=">"}
    $0 ~/^>/ {
        # match($2, /(.*).pdb*/)
        split($2,parts,"_")
        complex=""
        for (j = 1; j < length(parts); j++) {
            complex = complex parts[j]
            if (j < length(parts)-1){
                complex=complex"_" 
            }
        }
        if (!(complex in repComplex)) {
            print "#"complex
            repComplex[complex] = ""
        }
    }
    {print $0}
    ' "${1}" > "${1}.tmp" && mv -f -- "${1}.tmp" "${1}"
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
