#!/bin/sh -e

exists() {
	[ -f "$1" ]
}

if [ "${CLUSEARCH_PAR}" = 0 ]; then
    if exists "${IN}_ss.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" lndb "${IN}_h"  "${OUT}_tmp_ss_h" ${VERBOSITY} \
            || fail "lndb died"

        # shellcheck disable=SC2086
        "$MMSEQS" lndb "${IN}_ss"  "${OUT}_tmp_ss" ${VERBOSITY} \
            || fail "lndb died"

        # shellcheck disable=SC2086
        "$MMSEQS" base:makepaddedseqdb "${OUT}_tmp_ss" "${OUT}_ss" ${MAKEPADDEDSEQDB_PAR} \
            || fail "mmseqs makepaddedseqdb died"

        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${OUT}_tmp_ss" ${VERBOSITY} \
            || fail "rmdb died"

        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${OUT}_tmp_ss_h" ${VERBOSITY} \
            || fail "rmdb died"

        awk '{ print $3"\t"$1; }' "${OUT}_ss.lookup" > "${OUT}_ss.gpu_mapping1"


        if exists "${IN}.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" renamedbkeys "${OUT}_ss.gpu_mapping1" "${IN}" "${OUT}" \
                --subdb-mode 1 ${THREADS_PAR} \
                || fail "renamedbkeys died"
        fi

        if exists "${IN}_ca.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" renamedbkeys "${OUT}_ss.gpu_mapping1" "${IN}_ca" "${OUT}_ca" \
                --subdb-mode 1 ${THREADS_PAR} \
                || fail "renamedbkeys died"
        fi

        if exists "${IN}_h.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" renamedbkeys "${OUT}_ss.gpu_mapping1" "${IN}_h" "${OUT}_h" \
                --subdb-mode 1 ${THREADS_PAR} \
                || fail "renamedbkeys died"
        fi

        rm -f -- "${OUT}.lookup"
        awk '{print $1"\t"$2"\t"int($3/2)}' "${OUT}_ss.lookup" | sort -nk3 > "${OUT}.lookup"
        rm -f -- "${OUT}_ss.gpu_mapping1"

    else
        if exists "${IN}.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" base:makepaddedseqdb "${IN}" "${OUT}" ${MAKEPADDEDSEQDB_PAR} \
                || fail "mmseqs makepaddedseqdb died"
        fi

        if exists "${IN}_ca.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" base:makepaddedseqdb "${IN}_ca" "${OUT}_ca" ${MAKEPADDEDSEQDB_PAR} \
                || fail "mmseqs makepaddedseqdb died"
        fi
    fi
else
    if exists "${IN}_ss.dbtype"; then
        # # shellcheck disable=SC2086
        "$MMSEQS" lndb "${IN}_h"  "${OUT}_tmp_ss_h" ${VERBOSITY} \
            || fail "lndb died"

        # shellcheck disable=SC2086
        "$MMSEQS" lndb "${IN}_ss"  "${OUT}_tmp_ss" ${VERBOSITY} \
            || fail "lndb died"

        # shellcheck disable=SC2086
        "$MMSEQS" base:makepaddedseqdb "${OUT}_tmp_ss" "${OUT}_ss" ${MAKEPADDEDSEQDB_PAR} \
            || fail "mmseqs makepaddedseqdb died"

        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${OUT}_tmp_ss" ${VERBOSITY} \
            || fail "rmdb died"

        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${OUT}_tmp_ss_h" ${VERBOSITY} \
            || fail "rmdb died"

        awk '{ print $3"\t"$1; }' "${OUT}_ss.lookup" > "${OUT}_ss.gpu_mapping1"

        awk 'BEGIN{i=0} FNR==NR{name[$3]=1;print $3"\t"$1;i++; next} {if (!($1 in name)){print $1"\t"i; i++}}' \
        "${OUT}_ss.lookup" "${IN}.lookup" > "${OUT}_ss.gpu_mapping2"

        if exists "${IN}.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" renamedbkeys "${OUT}_ss.gpu_mapping1" "${IN}" "${OUT}" \
                --subdb-mode 1 ${THREADS_PAR} \
                || fail "renamedbkeys died"
        fi

        if exists "${IN}_ca.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" renamedbkeys "${OUT}_ss.gpu_mapping1" "${IN}_ca" "${OUT}_ca" \
                --subdb-mode 1 ${THREADS_PAR} \
                || fail "renamedbkeys died"
        fi

        if exists "${IN}_h.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" renamedbkeys "${OUT}_ss.gpu_mapping1" "${IN}_h" "${OUT}_h" \
                --subdb-mode 1 ${THREADS_PAR} \
                || fail "renamedbkeys died"
        fi

        # shellcheck disable=SC2086
        "$MMSEQS" renamedbkeys "${OUT}_ss.gpu_mapping2" "${IN}_seq" "${OUT}_seq" \
            --subdb-mode 1 ${THREADS_PAR} \
            || fail "renamedbkeys died"

        # shellcheck disable=SC2086
        "$MMSEQS" renamedbkeys "${OUT}_ss.gpu_mapping2" "${IN}_seq_ca" "${OUT}_seq_ca"  \
            --subdb-mode 1 ${THREADS_PAR} \
            || fail "renamedbkeys died"

        # shellcheck disable=SC2086
        "$MMSEQS" renamedbkeys "${OUT}_ss.gpu_mapping2" "${IN}_seq_h" "${OUT}_seq_h"  \
            --subdb-mode 1 ${THREADS_PAR} \
            || fail "renamedbkeys died"
        
        # shellcheck disable=SC2086
        "$MMSEQS" renamedbkeys "${OUT}_ss.gpu_mapping2" "${IN}_seq_ss" "${OUT}_seq_ss"  \
            --subdb-mode 1 ${THREADS_PAR} \
            || fail "renamedbkeys died"

        if exists "${IN}_clu.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" filterdb "${IN}_clu" "${OUT}_clutmp" --mapping-file  "${OUT}_ss.gpu_mapping2"  ${VERBOSITY} ${THREADS_PAR} \
                || fail "filterdb died"
            # shellcheck disable=SC2086
            "$MMSEQS" renamedbkeys "${OUT}_ss.gpu_mapping1" "${OUT}_clutmp" "${OUT}_clu"  \
                --subdb-mode 0 ${THREADS_PAR} \
                || fail "renamedbkeys died"  

            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${OUT}_clutmp" ${VERBOSITY} \
                || fail "rmdb died"
        fi

        if exists "${IN}_aln.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" filterdb "${IN}_aln" "${OUT}_alntmp" --mapping-file  "${OUT}_ss.gpu_mapping2"  ${VERBOSITY} ${THREADS_PAR} \
                || fail "filterdb died"
            # shellcheck disable=SC2086
            "$MMSEQS" renamedbkeys "${OUT}_ss.gpu_mapping1" "${OUT}_alntmp" "${OUT}_aln"  \
                --subdb-mode 0 ${THREADS_PAR} \
                || fail "renamedbkeys died"  
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${OUT}_alntmp" ${VERBOSITY} \
                || fail "rmdb died"
        fi

        rm -f -- "${OUT}.lookup"
        awk '{print $1"\t"$2"\t"int($3/2)}' "${OUT}_ss.lookup" | sort -nk3  > "${OUT}.lookup"
        rm -f -- "${OUT}_ss.gpu_mapping1"
        rm -f -- "${OUT}_ss.gpu_mapping2"

        rm "${OUT}_seq"
        rm "${OUT}_seq_ss"
        rm "${OUT}_seq_ca"
        rm "${OUT}_seq_h"

    else
        if exists "${IN}.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" base:makepaddedseqdb "${IN}" "${OUT}" ${MAKEPADDEDSEQDB_PAR} \
                || fail "mmseqs makepaddedseqdb died"
        fi
    fi
fi

if [ -e "${OUT}.sh" ]; then
  rm -f -- "${OUT}.sh"
fi 