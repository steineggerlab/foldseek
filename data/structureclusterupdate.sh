#!/bin/sh -e

fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
    [ ! -f "$1" ]
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

log() {
    if [ "${VERBOSITY}" = "-v 3" ]; then
        echo "$@"
    fi
}

# check number of input variables
[ "$#" -ne 6 ] && echo "Please provide <i:oldSequenceDB> <i:newSequenceDB> <i:oldClusteringDB> <o:newMappedSequenceDB> <o:newClusteringDB> <o:tmpDir>" && exit 1
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1
[ ! -f "$3.dbtype" ] && echo "$3.dbtype not found!" && exit 1
[   -f "$5.dbtype" ] && echo "$5.dbtype exists already!" && exit 1
[ ! -d "$6" ] && echo "tmp directory $6 not found!" && exit 1

OLDDB="$(abspath "$1")"
NEWDB="$(abspath "$2")"
OLDCLUST="$(abspath "$3")"
NEWMAPDB="$(abspath "$4")"
NEWCLUST="$(abspath "$5")"
TMP_PATH="$(abspath "$6")"

if notExists "${TMP_PATH}/removedSeqs"; then
    # shellcheck disable=SC2086
    "$MMSEQS" diffseqdbs "$OLDDB" "$NEWDB" "${TMP_PATH}/removedSeqs" "${TMP_PATH}/mappingSeqs" "${TMP_PATH}/newSeqs" ${DIFF_PAR} \
        || fail "Diff died"
fi

if [ ! -s "${TMP_PATH}/mappingSeqs" ]; then
    cat <<WARN
WARNING: There are no common sequences between $OLDDB and $NEWDB.
If you aim to add the sequences of $NEWDB to your previous clustering $OLDCLUST, you can run:

foldseek concatdbs "$OLDDB" "$NEWDB" "${OLDDB}.withNewSequences"
foldseek concatdbs "${OLDDB}_h" "${NEWDB}_h" "${OLDDB}.withNewSequences_h"
foldseek structureclusterupdate "$OLDDB" "${OLDDB}.withNewSequences" "$OLDCLUST" newMappedDB "$NEWCLUST" "${TMP_PATH}"
WARN
    rm -f "${TMP_PATH}/removedSeqs" "${TMP_PATH}/mappingSeqs" "${TMP_PATH}/newSeqs"
    exit 1
fi

if [ -s "${TMP_PATH}/removedSeqs" ]; then
    if [ -n "${RECOVER_DELETED}" ]; then
        log "=== Recover removed sequences"
        if notExists "${TMP_PATH}/OLDDB.removedMapping"; then
            HIGHESTID="$(awk '$1 > max { max = $1 } END { print max }' "${NEWDB}.index")"
            awk -v highest="$HIGHESTID" 'BEGIN { start=highest+1 } { printf("%s\t%.0f\n", $1, start); start=start+1; }' \
                "${TMP_PATH}/removedSeqs" > "${TMP_PATH}/OLDDB.removedMapping"
            cat "${TMP_PATH}/OLDDB.removedMapping" >> "${TMP_PATH}/mappingSeqs"
        fi

        if notExists "${TMP_PATH}/NEWDB.withOld.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" renamedbkeys "${TMP_PATH}/OLDDB.removedMapping" "${OLDDB}" "${TMP_PATH}/OLDDB.removedDb" --subdb-mode 1 ${VERBOSITY} \
                || fail "renamedbkeys died"
            # shellcheck disable=SC2086
            "$MMSEQS" concatdbs "$NEWDB" "${TMP_PATH}/OLDDB.removedDb" "${TMP_PATH}/NEWDB.withOld" --preserve-keys --threads 1 ${VERBOSITY} \
                || fail "concatdbs died"
            # shellcheck disable=SC2086
            "$MMSEQS" concatdbs "${NEWDB}_h" "${TMP_PATH}/OLDDB.removedDb_h" "${TMP_PATH}/NEWDB.withOld_h" --preserve-keys --threads 1 ${VERBOSITY} \
                || fail "concatdbs died"
        fi
        NEWDB="${TMP_PATH}/NEWDB.withOld"

        if [ -n "$REMOVE_TMP" ]; then
            echo "Remove temporary files 1/3"
            rm -f "${TMP_PATH}/OLDDB.removedMapping"
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/OLDDB.removedDb" ${VERBOSITY}
        fi
    else
        if notExists "${TMP_PATH}/REMOVEDMEMBERS.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" createsubdb "${TMP_PATH}/removedSeqs" "${OLDCLUST}" "${TMP_PATH}/REMOVEDMEMBERS" --subdb-mode 0 ${NOWARNINGS_PAR} \
                || fail "createsubdb died"
        fi

        if notExists "${TMP_PATH}/REMOVEDMEMBERS.withoutDeleted.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" filterdb "${TMP_PATH}/REMOVEDMEMBERS" "${TMP_PATH}/REMOVEDMEMBERS.withoutDeleted" --filter-file "${TMP_PATH}/removedSeqs" --positive-filter ${THREADS_PAR} \
                || fail "filterdb died"
        fi

        if notExists "${TMP_PATH}/REMOVEDMEMBERS.tsv"; then
            # shellcheck disable=SC2086
            "$MMSEQS" prefixid "${TMP_PATH}/REMOVEDMEMBERS.withoutDeleted" "${TMP_PATH}/REMOVEDMEMBERS.withoutDeleted.tsv" --tsv ${VERBOSITY} \
                || fail "prefixid died"
            awk '{ print $2; }' "${TMP_PATH}/REMOVEDMEMBERS.withoutDeleted.tsv" > "${TMP_PATH}/REMOVEDMEMBERS.tsv"
        fi

        if notExists "${TMP_PATH}/OLCLUST.withoutDeletedKeys.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" createsubdb "${TMP_PATH}/mappingSeqs" "${OLDCLUST}" "${TMP_PATH}/OLCLUST.withoutDeletedKeys" --subdb-mode 1 ${NOWARNINGS_PAR} \
                || fail "createsubdb died"
        fi

        if notExists "${TMP_PATH}/OLCLUST.withoutDeleted.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" filterdb "${TMP_PATH}/OLCLUST.withoutDeletedKeys" "${TMP_PATH}/OLCLUST.withoutDeleted" --filter-file "${TMP_PATH}/removedSeqs" --positive-filter ${THREADS_PAR} \
                || fail "filterdb died"
        fi
        OLDCLUST="${TMP_PATH}/OLCLUST.withoutDeleted"
    fi
fi

if notExists "${TMP_PATH}/newMappingSeqs"; then
    log "=== Update new sequences with old keys"
    MAXID="$(awk '$1 > max { max = $1 } END { print max }' "${OLDDB}.index" "${NEWDB}.index")"
    awk -v highest="$MAXID" 'BEGIN { start=highest+1 } { printf("%s\t%.0f\n", $1, start); start=start+1; }' \
        "${TMP_PATH}/newSeqs" > "${TMP_PATH}/newSeqs.mapped"
    awk '{ print $2"\t"$1 }' "${TMP_PATH}/mappingSeqs" > "${TMP_PATH}/mappingSeqs.reverse"
    cat "${TMP_PATH}/mappingSeqs.reverse" "${TMP_PATH}/newSeqs.mapped" > "${TMP_PATH}/newMappingSeqs"
    awk '{ print $2 }' "${TMP_PATH}/newSeqs.mapped" > "${TMP_PATH}/newSeqs"
fi

if notExists "${NEWMAPDB}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" renamedbkeys "${TMP_PATH}/newMappingSeqs" "${NEWDB}" "${NEWMAPDB}" ${VERBOSITY} \
        || fail "renamedbkeys died"
fi
if [ -f "${NEWDB}_ss.dbtype" ] && notExists "${NEWMAPDB}_ss.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" renamedbkeys "${TMP_PATH}/newMappingSeqs" "${NEWDB}_ss" "${NEWMAPDB}_ss" ${VERBOSITY} \
        || fail "renamedbkeys for 3Di died"
fi
if [ -f "${NEWDB}_h.dbtype" ] && notExists "${NEWMAPDB}_h.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" renamedbkeys "${TMP_PATH}/newMappingSeqs" "${NEWDB}_h" "${NEWMAPDB}_h" ${VERBOSITY} \
        || fail "renamedbkeys for headers died"
fi
NEWDB="${NEWMAPDB}"

NEWSEQ="${TMP_PATH}/newSeqs"
if [ -s "${TMP_PATH}/removedSeqs" ] && [ -z "${RECOVER_DELETED}" ]; then
    cat "${TMP_PATH}/REMOVEDMEMBERS.tsv" "${TMP_PATH}/newSeqs" > "${TMP_PATH}/newSeqs.withMembers"
    NEWSEQ="${TMP_PATH}/newSeqs.withMembers"
fi

if notExists "${TMP_PATH}/NEWDB.newSeqs.dbtype"; then
    log "=== Filter out new from old sequences"
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${NEWSEQ}" "$NEWDB" "${TMP_PATH}/NEWDB.newSeqs" ${VERBOSITY} --subdb-mode 1 \
        || fail "createsubdb died"
fi

if [ -f "${NEWDB}_ss.dbtype" ] && notExists "${TMP_PATH}/NEWDB.newSeqs_ss.dbtype"; then
    log "=== Filter out new from old 3Di sequences"
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${NEWSEQ}" "${NEWDB}_ss" "${TMP_PATH}/NEWDB.newSeqs_ss" ${VERBOSITY} --subdb-mode 1 \
        || fail "createsubdb for 3Di died"
fi

if notExists "${TMP_PATH}/OLDDB.repSeq.dbtype"; then
    log "=== Extract representative AA sequences"
    # Use createsubdb (not result2repseq): afdb50_clu members are keyed on afdb50_seq (241M),
    # but afdb50 only contains the 66M cluster reps. createsubdb uses the index keys directly.
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "$OLDCLUST" "$OLDDB" "${TMP_PATH}/OLDDB.repSeq" --subdb-mode 1 ${VERBOSITY} \
        || fail "createsubdb for rep AA sequences died"
fi

if [ -f "${OLDDB}_ss.dbtype" ] && notExists "${TMP_PATH}/OLDDB.repSeq_ss.dbtype"; then
    log "=== Extract representative 3Di sequences"
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "$OLDCLUST" "${OLDDB}_ss" "${TMP_PATH}/OLDDB.repSeq_ss" --subdb-mode 1 ${VERBOSITY} \
        || fail "createsubdb for rep 3Di sequences died"
fi

if notExists "${TMP_PATH}/newSeqsPref.dbtype"; then
    log "=== Prefilter new sequences against cluster representatives"
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefilter "${TMP_PATH}/NEWDB.newSeqs" "${TMP_PATH}/OLDDB.repSeq" \
        "${TMP_PATH}/newSeqsPref" ${PREFILTER_PAR} \
        || fail "Prefilter died"
fi

if notExists "${TMP_PATH}/newSeqsHits.dbtype"; then
    log "=== Structurally align new sequences against cluster representatives"
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" $ALIGNMENT_ALGO "${TMP_PATH}/NEWDB.newSeqs" "${TMP_PATH}/OLDDB.repSeq" \
        "${TMP_PATH}/newSeqsPref" "${TMP_PATH}/newSeqsHits" ${ALIGNMENT_PAR} \
        || fail "Structural alignment died"
fi

if notExists "${TMP_PATH}/newSeqsHits.best.dbtype"; then
    log "=== Keep best rep per new sequence (prevents one sequence joining multiple clusters)"
    # shellcheck disable=SC2086
    "$MMSEQS" filterdb "${TMP_PATH}/newSeqsHits" "${TMP_PATH}/newSeqsHits.best" \
        --extract-lines 1 ${THREADS_PAR} \
        || fail "filterdb best-hit died"
fi

if notExists "${TMP_PATH}/newSeqsHits.swapped.all.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" swapdb "${TMP_PATH}/newSeqsHits.best" "${TMP_PATH}/newSeqsHits.swapped.all" ${THREADS_PAR} \
        || fail "swapdb died"
    awk '$3 > 1 { print 1; exit; }' "${TMP_PATH}/newSeqsHits.swapped.all.index" > "${TMP_PATH}/newSeqsHits.swapped.hasHits"
fi

if [ -s "${TMP_PATH}/newSeqsHits.swapped.hasHits" ] && notExists "${TMP_PATH}/newSeqsHits.swapped.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" filterdb "${TMP_PATH}/newSeqsHits.swapped.all" "${TMP_PATH}/newSeqsHits.swapped" --trim-to-one-column ${THREADS_PAR} \
        || fail "filterdb died"
fi

UPDATEDCLUST="${TMP_PATH}/updatedClust"
if [ -f "${TMP_PATH}/newSeqsHits.swapped.dbtype" ]; then
    if notExists "${TMP_PATH}/updatedClust.dbtype"; then
        log "=== Merge structurally assigned sequences with previous clustering"
        # shellcheck disable=SC2086
        "$MMSEQS" mergedbs "$OLDCLUST" "${TMP_PATH}/updatedClust" "$OLDCLUST" "${TMP_PATH}/newSeqsHits.swapped" ${VERBOSITY} \
            || fail "mergedbs died"
    fi
else
    UPDATEDCLUST="$OLDCLUST"
fi

if notExists "${TMP_PATH}/toBeClusteredSeparately.dbtype"; then
    log "=== Extract unmapped sequences for independent clustering"
    awk '$3 == 1 {print $1}' "${TMP_PATH}/newSeqsHits.index" > "${TMP_PATH}/noHitSeqList"
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${TMP_PATH}/noHitSeqList" "$NEWDB" "${TMP_PATH}/toBeClusteredSeparately" ${VERBOSITY} --subdb-mode 1 \
        || fail "createsubdb of not hit seq. died"
fi

if notExists "${TMP_PATH}/newClusters.dbtype" && [ -s "${TMP_PATH}/toBeClusteredSeparately.index" ]; then
    log "=== Cluster unmapped sequences"
    # shellcheck disable=SC2086
    "$MMSEQS" cluster "${TMP_PATH}/toBeClusteredSeparately" "${TMP_PATH}/newClusters" "${TMP_PATH}/cluster" ${CLUST_PAR} \
        || fail "cluster of new seq. died"
fi

if [ -f "${TMP_PATH}/newClusters.dbtype" ]; then
    if notExists "$NEWCLUST"; then
        log "=== Merge updated clustering with new clusters"
        # shellcheck disable=SC2086
        "$MMSEQS" concatdbs "${UPDATEDCLUST}" "${TMP_PATH}/newClusters" "$NEWCLUST" --preserve-keys ${THREADS_PAR} \
            || fail "concatdbs died"
    fi
else
    # shellcheck disable=SC2086
    "$MMSEQS" mvdb "${UPDATEDCLUST}" "$NEWCLUST" ${VERBOSITY}
fi

if [ -n "$REMOVE_TMP" ]; then
    rm -f "${TMP_PATH}/newSeqs.mapped" "${TMP_PATH}/mappingSeqs.reverse" "${TMP_PATH}/newMappingSeqs"
    rm -f "${TMP_PATH}/noHitSeqList" "${TMP_PATH}/mappingSeqs" "${TMP_PATH}/newSeqs" "${TMP_PATH}/removedSeqs"
    rm -f "${TMP_PATH}/newSeqsHits.swapped.hasHits"

    if [ -n "${RECOVER_DELETED}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/NEWDB.withOld" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/NEWDB.withOld_h" ${VERBOSITY}
    else
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/OLCLUST.withoutDeletedKeys" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/OLCLUST.withoutDeleted" ${VERBOSITY}
    fi

    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/newSeqsHits.swapped" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/newClusters" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/newSeqsHits" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/newSeqsHits.best" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/newSeqsPref" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/toBeClusteredSeparately" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/NEWDB.newSeqs" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/newSeqsHits.swapped.all" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/OLDDB.repSeq" ${VERBOSITY}
    if [ -f "${TMP_PATH}/OLDDB.repSeq_ss.dbtype" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/OLDDB.repSeq_ss" ${VERBOSITY}
    fi
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/updatedClust" ${VERBOSITY}

    rm -rf "${TMP_PATH}/search" "${TMP_PATH}/cluster"
    rm -f "${TMP_PATH}/structureclusterupdate.sh"
fi
