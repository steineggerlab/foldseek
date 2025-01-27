#!/bin/sh -e
# Iterative sequence search workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

QUERYDB="$1"

STEP=0
# processing
if [ -z "$NUM_IT" ]; then
    NUM_IT=3
fi
while [ "$STEP" -lt "$NUM_IT" ]; do
    # call prefilter module
    if notExists "$TMP_PATH/pref_tmp_${STEP}.done"; then
        PARAM="PREFILTER_PAR_$STEP"
        eval TMP="\$$PARAM"
        TOOL="prefilter"
        if [ "$PREFMODE" = "UNGAPPED" ]; then
            TOOL="ungappedprefilter"
            PARAM="UNGAPPEDPREFILTER_PAR_$STEP"
            eval TMP="\$$PARAM"
        fi
        if [ $STEP -eq 0 ]; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" $TOOL "${QUERYDB}_ss" "${TARGET_PREFILTER}${INDEXEXT}" "$TMP_PATH/pref_${STEP}" ${TMP} \
                || fail "Prefilter died"
        else
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" $TOOL "${QUERYDB}_ss" "${TARGET_PREFILTER}${INDEXEXT}" "$TMP_PATH/pref_tmp_${STEP}" ${TMP} \
                || fail "Prefilter died"
        fi
        touch "$TMP_PATH/pref_tmp_${STEP}.done"
    fi

    if [ $STEP -ge 1 ]; then
        if notExists "$TMP_PATH/pref_$STEP.done"; then
            STEPONE=$((STEP-1))
            # shellcheck disable=SC2086
            "$MMSEQS" subtractdbs "$TMP_PATH/pref_tmp_${STEP}" "$TMP_PATH/aln_${STEPONE}" "$TMP_PATH/pref_${STEP}" $SUBSTRACT_PAR \
                || fail "Substract died"
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "$TMP_PATH/pref_tmp_${STEP}" ${VERBOSITY}
        fi
        touch "$TMP_PATH/pref_${STEP}.done"
    fi

	# call alignment module
	if notExists "$TMP_PATH/aln_tmp_${STEP}.done"; then
	    PARAM="ALIGNMENT_PAR_$STEP"
        eval TMP="\$$PARAM"

        if [ $STEP -eq 0 ]; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" "${ALIGNMENT_ALGO}" "${QUERYDB}" "${TARGET_ALIGNMENT}${INDEXEXT}" "$TMP_PATH/pref_${STEP}" "$TMP_PATH/aln_${STEP}" ${TMP} \
                || fail "Alignment died"
        else
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" "${ALIGNMENT_ALGO}" "${QUERYDB}" "${TARGET_ALIGNMENT}${INDEXEXT}" "$TMP_PATH/pref_${STEP}" "$TMP_PATH/aln_tmp_${STEP}" ${TMP} \
                || fail "Alignment died"
        fi
        touch "$TMP_PATH/aln_tmp_$STEP.done"
    fi

    if [ $STEP -gt 0 ]; then
        if notExists "$TMP_PATH/aln_$STEP.done"; then
            STEPONE=$((STEP-1))
            if [ $STEP -ne $((NUM_IT - 1)) ]; then
                # shellcheck disable=SC2086
                "$MMSEQS" mergedbs "${QUERYDB}" "$TMP_PATH/aln_${STEP}" "$TMP_PATH/aln_${STEPONE}" "$TMP_PATH/aln_tmp_${STEP}" ${VERBOSITY} \
                    || fail "Alignment died"
            else
                # shellcheck disable=SC2086
                "$MMSEQS" mergedbs "${QUERYDB}" "${RESULTS}" "$TMP_PATH/aln_${STEPONE}" "$TMP_PATH/aln_tmp_${STEP}" ${VERBOSITY} \
                    || fail "Alignment died"
            fi
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "$TMP_PATH/aln_${STEPONE}" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "$TMP_PATH/aln_tmp_${STEP}" ${VERBOSITY}
            touch "$TMP_PATH/aln_${STEP}.done"
        fi
    fi

    # create profiles
    if [ $STEP -ne $((NUM_IT - 1)) ]; then
        if notExists "$TMP_PATH/profile_${STEP}.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" result2profile "${QUERYDB}" "${TARGET_ALIGNMENT}${INDEXEXT}" "$TMP_PATH/aln_${STEP}" "$TMP_PATH/profile_${STEP}" ${PROFILE_PAR} \
                || fail "Create profile died"
        fi
    fi
	QUERYDB="$TMP_PATH/profile_${STEP}"
	STEP=$((STEP+1))
done

if [ -n "$REMOVE_TMP" ]; then
    STEP=0
    while [ "${STEP}" -lt "${NUM_IT}" ]; do
        if [ ${STEP} -gt 0 ]; then
            rm -f -- "$TMP_PATH/aln_${STEP}.done" "$TMP_PATH/pref_${STEP}.done"
        fi
        rm -f -- "$TMP_PATH/aln_tmp_${STEP}.done" "$TMP_PATH/pref_tmp_${STEP}.done"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_${STEP}" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/aln_${STEP}" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/profile_${STEP}" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/profile_${STEP}_ss" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/profile_${STEP}_h" ${VERBOSITY}
        STEP=$((STEP+1))
    done
    rm -f "$TMP_PATH/structureiterativesearch.sh"
fi

