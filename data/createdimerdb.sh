#!/bin/sh -e


#1. sharing monomers by changing lookup files to be possible to have redundant first column -> failed while expandmultimer & scoremultimer
#2. sharing monomers by changing index files to share second&third columns -> it reduces db memory , but it doesn't reduce computation time

# get list of each chain's id of multimers in chainidxlist (remove monomer chainids)
# get dimer's original chainid and name in ${OUT}.lookuptmp

if [ -e "${IN}.dbtype" ]; then
    awk 'FNR==NR{
        name[$1]=1; next
    } {
        if ($1 in name) {
            print $0
        }
    }' "${IN}.index" "${IN}.lookup" > "${TMP_PATH}/sublookup"

    awk -v firstout="${TMP_PATH}/lookuptmp" 'BEGIN {prev = -1; mon = 0; idx = -1; names = -1} {
        if(!($3 == prev)) {
            if (mon > 1) {
                split(idx, words, " ")
                split(names, namesplit, " ")
                for (i = 1; i <= mon - 1; i++) {
                    print words[i]
                    for (j = i + 1; j <= mon; j++) {
                        print words[i]"\t"namesplit[i] >> firstout
                        print words[j]"\t"namesplit[j] >> firstout
                    }
                }
                print words[mon]
            }
            prev = $3
            mon = 1
            idx = $1
            names = $2
        } else {
            mon ++
            idx = idx" "$1
            names = names" "$2
        }
    } END {
        if (mon > 1) {
            split(idx, words, " ")
            split(names, namesplit, " ")
            for (i = 1; i <= mon - 1; i++) {
                print words[i]
                for (j = i + 1; j <= mon; j++) {
                    print words[i]"\t"namesplit[i] >> firstout
                    print words[j]"\t"namesplit[j] >> firstout
                }
            }
            print words[mon]
        }
    }' "${TMP_PATH}/sublookup" > "${TMP_PATH}/chainidxlist"

    awk -v secondout="${TMP_PATH}/reallookup" '{
        print NR-1"\t"$2"\t"int((NR-1)/2) >> secondout
        print NR-1"\t"$1
    }' "${TMP_PATH}/lookuptmp" > "${TMP_PATH}/realIdx_originalIdx"

    awk '{print $3"\t"$2}' "${TMP_PATH}/reallookup" | awk -F'_' 'NR%2==1{OFS="_"; $NF=""; sub(/_$/, "", $0); print $0}' > "${TMP_PATH}/realsource" 

    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${TMP_PATH}/chainidxlist" "${IN}" "${TMP_PATH}/out" --subdb-mode 0 ${VERBOSITY_PAR} \
        || fail "createsubdb died"

    rm -- "${TMP_PATH}/out_h" "${TMP_PATH}/out_h.dbtype" "${TMP_PATH}/out_h.index"
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "${TMP_PATH}/chainidxlist" "${IN}_h" "${TMP_PATH}/out_h" --subdb-mode 0 ${VERBOSITY_PAR} \
        || fail "createsubdb died"

    rm -- "${TMP_PATH}/out.lookup" "${TMP_PATH}/out.source"
    mv "${TMP_PATH}/reallookup" "${OUT}.lookup"
    mv "${TMP_PATH}/realsource" "${OUT}.source"

    awk 'FNR==NR{name[$1]=$2"\t"$3; next}{print $1"\t"name[$2]}' "${TMP_PATH}/out_h.index" "${TMP_PATH}/realIdx_originalIdx" > "${OUT}_h.index"
    awk 'FNR==NR{name[$1]=$2"\t"$3; next}{print $1"\t"name[$2]}' "${TMP_PATH}/out.index" "${TMP_PATH}/realIdx_originalIdx" > "${OUT}.index"
    awk 'FNR==NR{name[$1]=$2"\t"$3; next}{print $1"\t"name[$2]}' "${TMP_PATH}/out_ss.index" "${TMP_PATH}/realIdx_originalIdx" > "${OUT}_ss.index"
    awk 'FNR==NR{name[$1]=$2"\t"$3; next}{print $1"\t"name[$2]}' "${TMP_PATH}/out_ca.index" "${TMP_PATH}/realIdx_originalIdx" > "${OUT}_ca.index"

    mv "${TMP_PATH}/out_h" "${OUT}_h"
    mv "${TMP_PATH}/out" "${OUT}"
    mv "${TMP_PATH}/out_ss" "${OUT}_ss"
    mv "${TMP_PATH}/out_ca" "${OUT}_ca"
    mv "${TMP_PATH}/out_h.dbtype" "${OUT}_h.dbtype"
    mv "${TMP_PATH}/out.dbtype" "${OUT}.dbtype"
    mv "${TMP_PATH}/out_ss.dbtype" "${OUT}_ss.dbtype"
    mv "${TMP_PATH}/out_ca.dbtype" "${OUT}_ca.dbtype"
fi

if [ -n "${REMOVE_TMP}" ]; then
    rm "${TMP_PATH}/chainidxlist"
    rm "${TMP_PATH}/lookuptmp"
    rm "${TMP_PATH}/realIdx_originalIdx"
    rm "${TMP_PATH}/out_ca.index"
    rm "${TMP_PATH}/out_h.index"
    rm "${TMP_PATH}/out.index"
    rm "${TMP_PATH}/out_ss.index"
fi
