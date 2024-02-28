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

# check number of input variables
[ "$#" -ne 3 ] && echo "Please provide <query> <out> <tmpDir>" && exit 1;
# check if files exist
[ ! -d "$3" ] && echo "tmp directory $3 not found!" && mkdir -p "$3";

# DOING : createdb
if notExists "${TMP_PATH}/input.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "${INPUT}" "${TMP_PATH}/input" ${CREATEDB_PAR} \
        || fail "input createdb died"
fi
INPUT="${TMP_PATH}/input"

# DOING : complexcluster 
if notExists "${TMP_PATH}/cmpl_db.dbtype"; then
    $MMSEQS complexcluster "${INPUT}" "${RESULT}" "${TMP_PATH}" "${COMPLEXCLUSTER_PAR}" \
        || fail "Complexcluster died"
fi
INPUTT="${TMP_PATH}/cmpl"

# DOING : make tsv file
if notExists "${TMP_PATH}/cluster.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${INPUTT}" "${INPUTT}" "${RESULT}" "${TMP_PATH}/cluster.tsv" ${THREADS_PAR} \
        || fail "Convert Alignments died"
fi

# FIXME : make rep_seq.fasta
# TODO: figure out how to represent complex sequences as a single fasta entry?
if notExists "${TMP_PATH}/rep_seq.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" result2repseq "${INPUTT}" "${RESULT}" "${TMP_PATH}/clu_rep" ${RESULT2REPSEQ_PAR} \
            || fail "Result2repseq  died"

    # shellcheck disable=SC2086
    "$MMSEQS" result2flat "${INPUTT}" "${INPUTT}" "${TMP_PATH}/clu_rep" "${TMP_PATH}/rep_seq.fasta" --use-fasta-header ${VERBOSITY_PAR} \
            || fail "result2flat died"
fi

# FIXME : make all_seq.fasta, and how ?
if notExists "${TMP_PATH}/all_seqs.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createseqfiledb "${INPUTT}" "${RESULT}" "${TMP_PATH}/clu_seqs" ${THREADS_PAR} \
            || fail "Result2repseq  died"

    # shellcheck disable=SC2086
    "$MMSEQS" result2flat "${INPUTT}" "${INPUTT}" "${TMP_PATH}/clu_seqs" "${TMP_PATH}/all_seqs.fasta" ${VERBOSITY_PAR} \
            || fail "result2flat died"
fi

mv "${TMP_PATH}/all_seqs.fasta"  "${RESULT}_all_seqs.fasta"
mv "${TMP_PATH}/rep_seq.fasta"  "${RESULT}_rep_seq.fasta"
mv "${TMP_PATH}/cluster.tsv"  "${RESULT}_cluster.tsv"

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
    "$MMSEQS" rmdb "${INPUTT}" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${RESULT}" ${VERBOSITY_PAR}
    rm -f "${TMP_PATH}/easycomplexcluster.sh"
fi