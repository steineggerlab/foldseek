#!/bin/sh -e

IN="$1"
OUT="$2"

outname=$(basename "${OUT}")
gpu_mapping1="${outname}_ss.gpu_mapping1"
gpu_mapping2="${outname}_ss.gpu_mapping2"

if [ ! -d "${TMP_PATH}" ]; then
    mkdir -p "${TMP_PATH}"
fi

exists() {
	[ -f "$1" ]
}

if [ "${CLUSEARCH_PAR}" = 0 ]; then
    if [ -e "${IN}_ss.dbtype" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" lndb "${IN}_h"  "${IN}_ss_h" ${VERBOSITY} \
            || fail "lndb died"

        # shellcheck disable=SC2086
        "$MMSEQS" base:makepaddedseqdb "${IN}_ss" "${OUT}_ss" ${MAKEPADDEDSEQDB_PAR} \
            || fail "mmseqs makepaddedseqdb died"

        awk '{ print $3"\t"$1; }' "${OUT}_ss.lookup" > "${TMP_PATH}/${gpu_mapping1}"

        if [ -e "${IN}.dbtype" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" renamedbkeys "${TMP_PATH}/${gpu_mapping1}" "${IN}" "${OUT}" \
                --subdb-mode 1 ${THREADS_PAR} \
                || fail "renamedbkeys died"
        fi

        if [ -e "${IN}_ca.dbtype" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" renamedbkeys "${TMP_PATH}/${gpu_mapping1}" "${IN}_ca" "${OUT}_ca" \
                --subdb-mode 1 ${THREADS_PAR} \
                || fail "renamedbkeys died"
        fi
    else
        if [ -e "${IN}.dbtype" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" base:makepaddedseqdb "${IN}" "${OUT}" ${MAKEPADDEDSEQDB_PAR} \
                || fail "mmseqs makepaddedseqdb died"
        fi

        if [ -e "${IN}_ca.dbtype" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" base:makepaddedseqdb "${IN}_ca" "${OUT}_ca" ${MAKEPADDEDSEQDB_PAR} \
                || fail "mmseqs makepaddedseqdb died"
        fi
    fi
else
    if [ -e "${IN}_ss.dbtype" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" lndb "${IN}_h" "${IN}_ss_h" ${VERBOSITY} \
            || fail "lndb died"

        # shellcheck disable=SC2086
        "$MMSEQS" base:makepaddedseqdb "${IN}_ss" "${OUT}_ss" ${MAKEPADDEDSEQDB_PAR} \
            || fail "mmseqs makepaddedseqdb died"

        awk '{ print $3"\t"$1; }' "${OUT}_ss.lookup" > "${TMP_PATH}/${gpu_mapping1}"

        awk 'BEGIN{i=0} FNR==NR{name[$3]=1;print $3"\t"$1;i++; next} {if (!($1 in name)){print $1"\t"i; i++}}' \
        "${OUT}_ss.lookup" "${IN}.lookup" > "${TMP_PATH}/${gpu_mapping2}"

        # shellcheck disable=SC2086
        "$MMSEQS" renamedbkeys "${TMP_PATH}/${gpu_mapping1}" "${IN}" "${OUT}" \
            --subdb-mode 1 ${THREADS_PAR} \
            || fail "renamedbkeys died"

        # shellcheck disable=SC2086
        "$MMSEQS" renamedbkeys "${TMP_PATH}/${gpu_mapping1}" "${IN}_ca" "${OUT}_ca" \
            --subdb-mode 1 ${THREADS_PAR} \
            || fail "renamedbkeys died"

        # shellcheck disable=SC2086
        "$MMSEQS" renamedbkeys "${TMP_PATH}/${gpu_mapping2}" "${IN}_seq" "${OUT}_seq" \
            --subdb-mode 1 ${THREADS_PAR} \
            || fail "renamedbkeys died"

        # shellcheck disable=SC2086
        "$MMSEQS" renamedbkeys "${TMP_PATH}/${gpu_mapping2}" "${IN}_seq_ca" "${OUT}_seq_ca"  \
            --subdb-mode 1 ${THREADS_PAR} \
            || fail "renamedbkeys died"

        # shellcheck disable=SC2086
        "$MMSEQS" renamedbkeys "${TMP_PATH}/${gpu_mapping2}" "${IN}_seq_h" "${OUT}_seq_h"  \
            --subdb-mode 1 ${THREADS_PAR} \
            || fail "renamedbkeys died"
        
        # shellcheck disable=SC2086
        "$MMSEQS" renamedbkeys "${TMP_PATH}/${gpu_mapping2}" "${IN}_seq_ss" "${OUT}_seq_ss"  \
            --subdb-mode 1 ${THREADS_PAR} \
            || fail "renamedbkeys died"

        if exists "${IN}_clu.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" filterdb "${IN}_clu" "${TMP_PATH}/${OUT}_clutmp" --mapping-file  "${TMP_PATH}/${gpu_mapping2}"  ${VERBOSITY} ${THREADS_PAR} \
                || fail "filterdb died"
            # shellcheck disable=SC2086
            "$MMSEQS" renamedbkeys "${TMP_PATH}/${gpu_mapping1}" "${TMP_PATH}/${OUT}_clutmp" "${OUT}_clu"  \
                --subdb-mode 0 ${THREADS_PAR} \
                || fail "renamedbkeys died"  
        fi

        if exists "${IN}_aln.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" filterdb "${IN}_aln" "${TMP_PATH}/${OUT}_alntmp" --mapping-file  "${TMP_PATH}/${gpu_mapping2}"  ${VERBOSITY} ${THREADS_PAR} \
                || fail "filterdb died"
            # shellcheck disable=SC2086
            "$MMSEQS" renamedbkeys "${TMP_PATH}/${gpu_mapping1}" "${TMP_PATH}/${OUT}_alntmp" "${OUT}_aln"  \
                --subdb-mode 0 ${THREADS_PAR} \
                || fail "renamedbkeys died"  

        fi

        rm "${OUT}_seq"
        rm "${OUT}_seq_ss"
        rm "${OUT}_seq_ca"
        rm "${OUT}_seq_h"

    else
        if [ -e "${IN}.dbtype" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" base:makepaddedseqdb "${IN}" "${OUT}" ${MAKEPADDEDSEQDB_PAR} \
                || fail "mmseqs makepaddedseqdb died"
        fi
    fi
fi

if [ -e "${IN}_ss.dbtype" ]; then
    rm -f -- "${TMP_PATH}/${gpu_mapping1}"
    rm -f -- "${TMP_PATH}/${gpu_mapping2}"
    
    if [ -e "${IN}_aln.dbtype" ]; then
        "$MMSEQS" rmdb "${TMP_PATH}/${OUT}_alntmp" 
    fi
    if [ -e "${IN}_clu.dbtype" ]; then
        "$MMSEQS" rmdb "${TMP_PATH}/${OUT}_clutmp" 
    fi
fi

rm -f "${TMP_PATH}/makepaddeddb.sh"