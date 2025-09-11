#!/bin/sh -e


exists() {
	[ -f "$1" ]
}

if [ -e "${IN}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" filterdimerdb "${IN}" "${TMP_PATH}/contactlist" ${FILTERDIMERDB_PAR} \
        || fail "filterdimerdb died"
    sort -nk2 "${TMP_PATH}/contactlist.index" > "${TMP_PATH}/contactlist_2"
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "${TMP_PATH}/contactlist_2" "${IN}" "${OUT}" --subdb-mode 0 ${VERBOSITY_PAR} \
        || fail "createsubdb died"
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "${TMP_PATH}/contactlist_2" "${IN}_ss" "${OUT}_ss" --subdb-mode 0 ${VERBOSITY_PAR} \
        || fail "createsubdb died"
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "${TMP_PATH}/contactlist_2" "${IN}_ca" "${OUT}_ca" --subdb-mode 0 ${VERBOSITY_PAR} \
        || fail "createsubdb died"
    if exists "${IN}_id"; then
        # shellcheck disable=SC2086
        "$MMSEQS" base:createsubdb "${TMP_PATH}/contactlist_2" "${IN}_id" "${OUT}_id" --subdb-mode 0 ${VERBOSITY_PAR} \
            || fail "createsubdb died"
    fi
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${OUT}_h"
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${TMP_PATH}/contactlist_2" "${IN}_h" "${OUT}_h" \
        || fail "createsubdb died"
    
    sort -nk2 "${OUT}.index" | awk '{print NR-1"\t"$2"\t"$3}' > "${TMP_PATH}/db.index2"
    sort -nk2 "${OUT}_ss.index" | awk '{print NR-1"\t"$2"\t"$3}' > "${TMP_PATH}/db_ss.index2"
    sort -nk2 "${OUT}_ca.index" | awk '{print NR-1"\t"$2"\t"$3}' > "${TMP_PATH}/db_ca.index2"
    sort -nk2 "${OUT}_h.index" | awk '{print NR-1"\t"$2"\t"$3}' > "${TMP_PATH}/db_h.index2"
    sort -nk2 "${OUT}_id.index" | awk '{print NR-1"\t"$2"\t"$3}' > "${TMP_PATH}/db_id.index2"

    mv "${TMP_PATH}/db.index2" "${OUT}.index"
    mv "${TMP_PATH}/db_ss.index2" "${OUT}_ss.index"
    mv "${TMP_PATH}/db_ca.index2" "${OUT}_ca.index"
    mv "${TMP_PATH}/db_h.index2" "${OUT}_h.index"
    mv "${TMP_PATH}/db_id.index2" "${OUT}_id.index"

    rm "${OUT}.lookup" "${OUT}.source"
    awk 'FNR==NR{name[$1]=$2; source[$1]=$3; next} BEGIN{num=0}{print num"\tDI"int(num/2)"_"name[$1]"\t"int(num/2)"\t"source[$1]; num++}' "${IN}.lookup" "${TMP_PATH}/contactlist_2" > "${TMP_PATH}/lookuptmp"
    awk 'FNR==NR{name[$1]=$2; next} NR%2==1{print $3"\tDI"$3"_"name[$4]}' "${IN}.source" "${TMP_PATH}/lookuptmp" > "${OUT}.source"
    cut -f1,2,3 "${TMP_PATH}/lookuptmp" > "${OUT}.lookup"
    
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${OUT}" "${OUT}_h" "${TMP_PATH}/tmpheader" --threads 1 ${VERBOSITY_PAR}
    paste "${OUT}.lookup" "${TMP_PATH}/tmpheader" | awk -F"\t" '{print $1"\t"$2"\t"$5}' > "${TMP_PATH}/tmpheader2"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${OUT}_h"
    # shellcheck disable=SC2086
    "$MMSEQS" tsv2db "${TMP_PATH}/tmpheader2" "${OUT}_h" ${VERBOSITY_PAR}

    if exists "${OUT}_mapping"; then
        rm "${OUT}_mapping"
        awk 'FNR==NR{map[$1]=$2; next}BEGIN{num=0}{print num"\t"map[$1]; num++}' "${IN}_mapping" "${TMP_PATH}/contactlist_2"  > "${OUT}_mapping"
    fi

    if exists "${OUT}_taxonomy"; then
        rm "${OUT}_taxonomy"
        #TODO: taxonomy
    fi
fi

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/contactlist"
    rm "${TMP_PATH}/contactlist_2"
    rm "${TMP_PATH}/lookuptmp"
    rm "${TMP_PATH}/tmpheader"
    rm "${TMP_PATH}/tmpheader2"
    rm -f "${TMP_PATH}/createdimerdb.sh"
fi
