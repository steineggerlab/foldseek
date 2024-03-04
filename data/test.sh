#!/bin/sh -e


mapCmpl2Chain() {
    awk -F"\t" 'BEGIN {}
        NR==FNR { 
            if ($0 ~ /^[0-9]+$/) {
                reps[$1]=$1;next
            }
        }
        { if ($3 in reps) {
            print $1
            print "\0"$1
            }
        }
  ' "${1}" "${2}".lookup > "${3}"
}

# if notExists "${TMP_PATH}/complex_rep_seq.fasta"; then
#     mapCmpl2Chain "${TMP_PATH}/complex_clust" "${INPUT}" "${TMP_PATH}/complex_clust_chains"
#     # shellcheck disable=SC2086
#     "$MMSEQS" result2repseq "${INPUT}" "${TMP_PATH}/complex_clust_chains" "${TMP_PATH}/complex_clu_rep" ${RESULT2REPSEQ_PAR} \
#             || fail "Result2repseq  died"

#     # shellcheck disable=SC2086
#     "$MMSEQS" result2flat "${INPUT}" "${INPUT}"  "${TMP_PATH}/complex_clust_rep" "${TMP_PATH}/complex_rep_seq.fasta" --use-fasta-header ${VERBOSITY_PAR} \
#             || fail "result2flat died"
# fi

lookupfile="/Users/steineggerlab/Desktop/foldseek/toydata/emmanuel_db/emmanuel_7soy_db"
testfile="/Users/steineggerlab/Desktop/foldseek/toydata/emmanuel_7soy_easycc_tmp/9595115010967635338/complex_clust"
mapCmpl2Chain $testfile $lookupfile "/Users/steineggerlab/Desktop/foldseek/toydata/emmanuel_7soy_easycc_tmp/9595115010967635338/complex_clust_reps"
