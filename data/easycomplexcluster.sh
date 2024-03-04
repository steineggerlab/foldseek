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
    "$MMSEQS" complexcluster "${TMP_PATH}/input" "${TMP_PATH}/complex_clu" "$(dirname "${RESULT}")" ${COMPLEXCLUSTER_PAR} \
        || fail "Complexcluster died"
fi

INPUT="$(dirname "${RESULT}")/latest/complex_db"
if notExists "${TMP_PATH}/cluster.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${INPUT}" "${INPUT}" "${TMP_PATH}/complex_clu" "${TMP_PATH}/cluster.tsv" ${THREADS_PAR} \
        || fail "Convert Alignments died"
fi

#TODO: figure out how to represent complex sequences as a single fasta entry?
if notExists "${TMP_PATH}/complex_rep_seq.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" result2repseq "${INPUT}" "${TMP_PATH}/complex_clu" "${TMP_PATH}/complex_clu_rep" ${RESULT2REPSEQ_PAR} \
            || fail "Result2repseq  died"

    # shellcheck disable=SC2086
    "$MMSEQS" result2flat "${INPUT}" "${INPUT}"  "${TMP_PATH}/complex_clu_rep" "${TMP_PATH}/complex_rep_seq.fasta" --use-fasta-header ${VERBOSITY_PAR} \
            || fail "result2flat died"
fi

if notExists "${TMP_PATH}/complex_all_seqs.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createseqfiledb "${INPUT}" "${TMP_PATH}/complex_clu" "${TMP_PATH}/complex_clu_seqs" ${THREADS_PAR} \
            || fail "Result2repseq  died"

    # shellcheck disable=SC2086
    "$MMSEQS" result2flat "${INPUT}" "${INPUT}" "${TMP_PATH}/complex_clu_seqs" "${TMP_PATH}/complex_all_seqs.fasta" ${VERBOSITY_PAR} \
            || fail "result2flat died"
fi

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
    rm -f "${TMP_PATH}/easycomplexcluster.sh"
fi