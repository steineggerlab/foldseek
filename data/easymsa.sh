#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# Generate fake prefilter file for all vs all alignments
fake_pref() {
    QDB="$1"
    TDB="$2"
    RES="$3"
    # create link to data file which contains a list of all targets that should be aligned
    # ln -s "${TDB}.index" "${RES}"
	# TODO structurealign complains about no data file being found when using symlink
	cp "${TDB}.index" "${RES}"
    # create new index repeatedly pointing to same entry
    # INDEX_SIZE="$(echo $(wc -c < "${TDB}.index"))"
    INDEX_SIZE=$(wc -c < "${TDB}.index")
    awk -v size="${INDEX_SIZE}" '{ print $1"\t0\t"size; }' "${QDB}.index" > "${RES}.index"
    # create dbtype (7)
    awk 'BEGIN { printf("%c%c%c%c",7,0,0,0); exit; }' > "${RES}.dbtype"
}

STRUCTUREDB="${TMP_PATH}/structures"
if notExists "${STRUCTUREDB}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "$@" "${STRUCTUREDB}" ${CREATEDB_PAR} \
        || fail "Structure createdb died"
fi

if notExists "${TMP_PATH}/perf.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" kmermatcher "${STRUCTUREDB}_ss" "${TMP_PATH}/perf" ${KMERMATCHER_PAR} \
        || fail "Structure kmermatcher died"
fi

if notExists "${TMP_PATH}/aln.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" structurealign "${STRUCTUREDB}" "${STRUCTUREDB}" "${TMP_PATH}/perf" "${TMP_PATH}/aln" ${ALIGN_PAR} \
        || fail "Structure align died"
fi

if notExists "${TMP_PATH}/clu.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" clust "${STRUCTUREDB}" "${TMP_PATH}/aln" "${TMP_PATH}/clu" ${CLUST_PAR} \
        || fail "Structure clust died"
fi

if notExists "${TMP_PATH}/${TREE}"; then
    # shellcheck disable=SC2086
	# Query DB, Alignment DB, temporary directory
	"$MMSEQS" structuremsa "${STRUCTUREDB}" "${TMP_PATH}/clu" "${RESULTS}" ${STRUCTUREMSA_PAR} \
		|| fail "structuremsa died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    if [ -n "${GREEDY_BEST_HITS}" ]; then
        # shellcheck disable=SC2086
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/result_best" ${VERBOSITY}
    fi
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result" ${VERBOSITY}
    if [ -z "${LEAVE_INPUT}" ]; then
        if [ -f "${TMP_PATH}/target" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_h" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_ca" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_ss" ${VERBOSITY}
        fi
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_h" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_ca" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_ss" ${VERBOSITY}
    fi
    rm -rf "${TMP_PATH}/search_tmp"
    rm -f "${TMP_PATH}/easymsa.sh"
fi
