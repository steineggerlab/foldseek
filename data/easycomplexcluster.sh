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


if notExists "${TMP_PATH}/input.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "${INPUT}" "${TMP_PATH}/input" ${CREATEDB_PAR} \
        || fail "input createdb died"
fi

if notExists "${TMP_PATH}/complex_clu.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" complexcluster "${TMP_PATH}/input" "${TMP_PATH}/complex_clu" "${TMP_PATH}" ${COMPLEXCLUSTER_PAR} \
        || fail "Complexcluster died"
fi

INPUT="${TMP_PATH}/latest/complex_db"
if notExists "${TMP_PATH}/cluster.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${INPUT}" "${INPUT}" "${TMP_PATH}/complex_clust" "${TMP_PATH}/cluster.tsv" ${THREADS_PAR} \
        || fail "Convert Alignments died"
fi

#TODO: move it to complexcluster.sh?

mapCmpl2Chain() {
    awk 'BEGIN {FS="\t"}
        NR==FNR {
            if (!($0 ~ /^[0-9]+$/)) {
                split($1,name,"\0")
                reps[name[2]]=$1;next
            } else if (FNR==1) {
                reps[$1]=$1
            }
            next
        }
        { if ($3 in reps) {
            print "\0"$1
        }
    }' "${1}" "${2}".lookup > "${3}"
}

if notExists "${TMP_PATH}/complex_rep_seq.fasta"; then
    mapCmpl2Chain "${TMP_PATH}/complex_clust" "${INPUT}" "${TMP_PATH}/complex_clust_chains"
    # shellcheck disable=SC2086
    "$MMSEQS" result2repseq "${INPUT}" "${TMP_PATH}/complex_clust_chains" "${TMP_PATH}/complex_clust_rep" ${RESULT2REPSEQ_PAR} \
            || fail "Result2repseq  died"

    # shellcheck disable=SC2086
    "$MMSEQS" result2flat "${INPUT}" "${INPUT}"  "${TMP_PATH}/complex_clust_rep" "${TMP_PATH}/complex_rep_seq.fasta" --use-fasta-header ${VERBOSITY_PAR} \
            || fail "result2flat died"
fi

# if notExists "${TMP_PATH}/complex_all_seqs.fasta"; then
#     # shellcheck disable=SC2086
#     "$MMSEQS" createseqfiledb "${INPUT}" "${TMP_PATH}/complex_clust" "${TMP_PATH}/complex_clust_seqs" ${THREADS_PAR} \
#             || fail "Result2repseq  died"

#     # shellcheck disable=SC2086
#     "$MMSEQS" result2flat "${INPUT}" "${INPUT}" "${TMP_PATH}/complex_clust_seqs" "${TMP_PATH}/complex_all_seqs.fasta" ${VERBOSITY_PAR} \
#             || fail "result2flat died"
# fi

mv "${TMP_PATH}/complex_all_seqs.fasta"  "${RESULT}_all_seqs.fasta"
mv "${TMP_PATH}/complex_rep_seq.fasta"  "${RESULT}_rep_seq.fasta"
mv "${TMP_PATH}/cluster.tsv"  "${RESULT}_cluster.tsv"

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input_h" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/complex_db" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/complex_clu_seqs" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/complex_clu_rep" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/complex_clu" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${INPUT}" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${INPUT}_h" ${VERBOSITY_PAR}
    rm -rf "${TMP_PATH}/latest"
    rm -f "${TMP_PATH}/easycomplexcluster.sh"
fi