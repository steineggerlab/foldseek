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

# check if files exist
[ ! -f "${INPUT}.dbtype" ] && echo "${INPUT}.dbtype not found!" && exit 1;
[ ! -d "${TMP_PATH}" ] && echo "tmp directory ${TMP_PATH} not found!" && mkdir -p "${TMP_PATH}";

# DOING : createdb
if notExists "${INPUT}.dbtype"; then
    if notExists "${TMP_PATH}/query"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createdb "${INPUT}" "${TMP_PATH}/input" ${CREATEDB_PAR} \
            || fail "input createdb died"
    fi
fi

# DOING : complexsearch
if notExists "${TMP_PATH}/complex_result.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" complexsearch "${INPUT}" "${INPUT}" "${TMP_PATH}/complex_result" "${TMP_PATH}/complexsearch_tmp" ${COMPLEXSEARCH_PAR} \
        || fail "ComplexSearch died"
fi
COMPDB="${TMP_PATH}/complexsearch_tmp"

# DOING : call complexcluster or filtercomplex+awk
# TODO : maybe save filtercomplex result file to sub-dir of TMP_PATH
if notExists "${TMP_PATH}/${RESULT}.dbtype"; then
    $MMSEQS "${CLUSTER_MODULE}" "${INPUT}" "${INPUT}" "${TMP_PATH}/${RESULT}" "${TMP_PATH}" ${COMPLEXCLUSTER_PAR} \
        || fail "ClusterComplex died"
fi

# DOING : make tsv file
if notExists "${TMP_PATH}/cluster.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${INPUT}" "${INPUT}" "${TMP_PATH}/${RESULTS}" "${TMP_PATH}/cluster.tsv" ${THREADS_PAR} \
        || fail "Convert Alignments died"
fi

# FIXME : make rep_seq.fasta, and how ?
# TODO: figure out how to represent complex sequences as a single fasta entry?
if notExists "${TMP_PATH}/rep_seq.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" result2repseq "${INPUT}" "${TMP_PATH}/${RESULTS}" "${TMP_PATH}/clu_rep" ${RESULT2REPSEQ_PAR} \
            || fail "Result2repseq  died"

    # shellcheck disable=SC2086
    "$MMSEQS" result2flat "${INPUT}" "${INPUT}" "${TMP_PATH}/clu_rep" "${TMP_PATH}/rep_seq.fasta" --use-fasta-header ${VERBOSITY_PAR} \
            || fail "result2flat died"
fi

# FIXME : make all_seq.fasta, and how ?
if notExists "${TMP_PATH}/all_seqs.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createseqfiledb "${INPUT}" "${TMP_PATH}/${RESULTS}" "${TMP_PATH}/clu_seqs" ${THREADS_PAR} \
            || fail "Result2repseq  died"

    # shellcheck disable=SC2086
    "$MMSEQS" result2flat "${INPUT}" "${INPUT}" "${TMP_PATH}/clu_seqs" "${TMP_PATH}/all_seqs.fasta" ${VERBOSITY_PAR} \
            || fail "result2flat died"
fi

mv "${TMP_PATH}/all_seqs.fasta"  "${RESULTS}_all_seqs.fasta"
mv "${TMP_PATH}/rep_seq.fasta"  "${RESULTS}_rep_seq.fasta"
mv "${TMP_PATH}/cluster.tsv"  "${RESULTS}_cluster.tsv"

# TODO : remove tmp -> tide up and organize
if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input_h" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/clu_seqs" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/clu_rep" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$2" ${VERBOSITY_PAR}
    rm -rf "${TMP_PATH}/complexsearch_tmp"
    rm -f "${TMP_PATH}/easycomplexcluster.sh"
fi