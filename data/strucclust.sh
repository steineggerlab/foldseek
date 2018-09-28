#!/bin/sh -e
# Assembler workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# check input variables
[ ! -n "${OUT_FILE}" ] && echo "Please provide OUT_FILE" && exit 1
[ ! -n "${TMP_PATH}" ] && echo "Please provide TMP_PATH" && exit 1

# check if files exists
[   -f "${OUT_FILE}" ] &&  echo "${OUT_FILE} exists already!" && exit 1
[ ! -d "${TMP_PATH}" ] &&  echo "tmp directory ${TMP_PATH} not found!" && mkdir -p "${TMP_PATH}"

if notExists "${TMP_PATH}/nucl_reads"; then
    if [ -n "${PAIRED_END}" ]; then
        echo "PAIRED END MODE"
        # shellcheck disable=SC2086
        "$MMSEQS" mergereads "$@" "${TMP_PATH}/nucl_reads" ${VERBOSITY_PAR} \
            || fail "mergereads failed"
    else
        # shellcheck disable=SC2086
        "$MMSEQS" createdb "$@" "${TMP_PATH}/nucl_reads" ${CREATEDB_PAR} \
            || fail "createdb failed"
    fi
fi

INPUT="${TMP_PATH}/nucl_reads"
if notExists "${TMP_PATH}/nucl_6f_start"; then
    # shellcheck disable=SC2086
    "$MMSEQS" extractorfs "${INPUT}" "${TMP_PATH}/nucl_6f_start" ${EXTRACTORFS_START_PAR} \
        || fail "extractorfs start step died"
fi

if notExists "${TMP_PATH}/aa_6f_start"; then
    # shellcheck disable=SC2086
    "$MMSEQS" translatenucs "${TMP_PATH}/nucl_6f_start" "${TMP_PATH}/aa_6f_start" ${TRANSLATENUCS_PAR} \
        || fail "translatenucs start step died"
fi

if notExists "${TMP_PATH}/nucl_6f_long"; then
    # shellcheck disable=SC2086
    "$MMSEQS" extractorfs "${INPUT}" "${TMP_PATH}/nucl_6f_long" ${EXTRACTORFS_LONG_PAR} \
        || fail "extractorfs longest step died"
fi

if notExists "${TMP_PATH}/aa_6f_long"; then
    # shellcheck disable=SC2086
    "$MMSEQS" translatenucs "${TMP_PATH}/nucl_6f_long" "${TMP_PATH}/aa_6f_long" ${TRANSLATENUCS_PAR} \
        || fail "translatenucs long step died"
fi

if notExists "${TMP_PATH}/aa_6f_start_long"; then
    # shellcheck disable=SC2086
    "$MMSEQS" concatdbs "${TMP_PATH}/aa_6f_long" "${TMP_PATH}/aa_6f_start" "${TMP_PATH}/aa_6f_start_long" ${VERBOSITY_PAR} \
        || fail "concatdbs start long step died"
fi

if notExists "${TMP_PATH}/aa_6f_start_long_h"; then
    # shellcheck disable=SC2086
    "$MMSEQS" concatdbs "${TMP_PATH}/nucl_6f_long_h" "${TMP_PATH}/nucl_6f_start_h" "${TMP_PATH}/aa_6f_start_long_h" ${VERBOSITY_PAR} \
        || fail "concatdbs start long step died"
fi

INPUT="${TMP_PATH}/aa_6f_start_long"
STEP=0
if [ -z "$NUM_IT" ]; then
    NUM_IT=1
fi

while [ "$STEP" -lt "$NUM_IT" ]; do
    echo "STEP: $STEP"

    # 1. Finding exact $k$-mer matches.
    if notExists "${TMP_PATH}/pref_$STEP"; then
        PARAM=KMERMATCHER${STEP}_PAR
        eval KMERMATCHER_TMP="\$$PARAM"
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" kmermatcher "$INPUT" "${TMP_PATH}/pref_$STEP" ${KMERMATCHER_TMP} \
            || fail "Kmer matching step died"
    fi

    # 2. Ungapped alignment
    if notExists "${TMP_PATH}/aln_$STEP"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" rescorediagonal "$INPUT" "$INPUT" "${TMP_PATH}/pref_$STEP" "${TMP_PATH}/aln_$STEP" ${UNGAPPED_ALN_PAR} \
            || fail "Ungapped alignment step died"
    fi

    ALN="${TMP_PATH}/aln_$STEP"
    if [ $STEP -eq 0 ]; then
        if notExists "${TMP_PATH}/corrected_seqs"; then
            # shellcheck disable=SC2086
            "$MMSEQS" findassemblystart "$INPUT" "${TMP_PATH}/aln_$STEP" "${TMP_PATH}/corrected_seqs" ${THREADS_PAR} \
                || fail "Findassemblystart alignment step died"
        fi
        INPUT="${TMP_PATH}/corrected_seqs"
        if notExists "${TMP_PATH}/aln_corrected_$STEP"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" rescorediagonal "$INPUT" "$INPUT" "${TMP_PATH}/pref_$STEP" "${TMP_PATH}/aln_corrected_$STEP" ${UNGAPPED_ALN_PAR} \
                || fail "Ungapped alignment step died"
        fi
        ALN="${TMP_PATH}/aln_corrected_$STEP"
    fi

    # 3. Assemble
    if notExists "${TMP_PATH}/assembly_$STEP"; then
        # shellcheck disable=SC2086
        "$MMSEQS" assembleresults "$INPUT" "${ALN}" "${TMP_PATH}/assembly_$STEP" ${ASSEMBLE_RESULT_PAR} \
            || fail "Assembly step died"
    fi

    INPUT="${TMP_PATH}/assembly_$STEP"
    STEP="$((STEP+1))"
done
STEP="$((STEP-1))"

# post processing
RESULT="${TMP_PATH}/assembly_${STEP}"
if [ -n "${PROTEIN_FILTER}" ]; then
    RESULT="${TMP_PATH}/assembly_${STEP}_filtered"
    if notExists "${TMP_PATH}/assembly_${STEP}_filtered"; then
        # shellcheck disable=SC2086
        "$MMSEQS" filternoncoding "${TMP_PATH}/assembly_${STEP}" "${TMP_PATH}/assembly_${STEP}_filtered" ${FILTERNONCODING_PAR} \
            || fail "filternoncoding died"
    fi
fi

# select only assembled sequences
if notExists "${RESULT}_only_assembled.index"; then
    # detect assembled proteins sequences
    awk 'NR == FNR { f[$1] = $0; next } $1 in f { print f[$1], $0 }' "${RESULT}.index" "${TMP_PATH}/aa_6f_start_long.index" > "${RESULT}_tmp.index"
    awk '$3 > $6 { print $1"\t"$2"\t"$3 }' "${RESULT}_tmp.index" > "${RESULT}_only_assembled1.index"
    # detect complete proteins with * at start and end
    awk '/^\x00?\*[A-Z]*\*$/{ f[NR-1]=1; next } $1 in f { print $0 }' "${RESULT}" "${RESULT}.index" > "${RESULT}_only_assembled2.index"
    # keep only non-redundant entries
    cat "${RESULT}_only_assembled1.index" "${RESULT}_only_assembled2.index" | sort | uniq > "${RESULT}_only_assembled.index"
fi

# create fasta output
if notExists "${RESULT}_only_assembled"; then
    ln -s "${RESULT}" "${RESULT}_only_assembled"
fi

if notExists "${RESULT}_only_assembled_h"; then
    ln -s "${TMP_PATH}/aa_6f_start_long_h" "${RESULT}_only_assembled_h"
fi

if notExists "${RESULT}_only_assembled_h.index"; then
    ln -s "${TMP_PATH}/aa_6f_start_long_h.index" "${RESULT}_only_assembled_h.index"
fi

if notExists "${RESULT}_only_assembled.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" convert2fasta "${RESULT}_only_assembled" "${RESULT}_only_assembled.fasta" ${VERBOSITY_PAR} \
        || fail "convert2fasta died"
fi

mv -f "${RESULT}_only_assembled.fasta" "$OUT_FILE" \
    || fail "Could not move result to $OUT_FILE"

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    rm -f "${TMP_PATH}/pref_"*
    rm -f "${TMP_PATH}/aln_"*
    rm -f "${TMP_PATH}/assembly_"*
fi
