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

while [ "$STEP" -lt "$NUM_IT" ]; do
  if notExists "$TMP_PATH/pref_${STEP}.dbtype"; then
      PARAM="PREFILTER_PAR_$STEP"
      eval TMP="\$$PARAM"
      TOOL="prefilter"
      if [ "$PREFMODE" = "UNGAPPED" ]; then
          TOOL="ungappedprefilter"
          PARAM="UNGAPPEDPREFILTER_PAR_$STEP"
          eval TMP="\$$PARAM"
      fi
      # shellcheck disable=SC2086
      $RUNNER "$MMSEQS" $TOOL "${QUERYDB}_ss" "${TARGET_PREFILTER}${INDEXEXT}" "$TMP_PATH/pref_${STEP}" ${TMP} \
          || fail "Prefilter died"
  fi

  # call alignment module
  PARAM="ALIGNMENT_PAR_$STEP"
  eval TMP="\$$PARAM"
  if [ $STEP -eq 0 ]; then
    if notExists "$TMP_PATH/aln_${STEP}.dbtype"; then
      # shellcheck disable=SC2086
      $RUNNER "$MMSEQS" "${ALIGNMENT_ALGO}" "${QUERYDB}" "${TARGET_ALIGNMENT}${INDEXEXT}" "$TMP_PATH/pref_${STEP}" "$TMP_PATH/aln_${STEP}" ${TMP} \
            || fail "Alignment died"
    fi
  else
    INTERMEDIATE="$TMP_PATH/pref_${STEP}"
    if [ -n "${EXPAND}" ]; then
        if notExists "${TMP_PATH}/aln_${STEP}.dbtype"; then
          # shellcheck disable=SC2086
          $RUNNER "$MMSEQS" "${ALIGNMENT_ALGO}" "${QUERYDB}" "${TARGET_ALIGNMENT}${INDEXEXT}" "$TMP_PATH/pref_${STEP}" "${TMP_PATH}/aln_${STEP}" ${TMP} \
                   || fail "Alignment died"
        fi
        if notExists "${TMP_PATH}/aln_expanded.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" mergeresultsbyset "${TMP_PATH}/aln_${STEP}" "${TARGET_ALIGNMENT}${INDEXEXT}" "${TMP_PATH}/aln_expanded" ${MERGERESULTBYSET_PAR} \
                || fail "Expand died"
            # shellcheck disable=SC2086
            "$MMSEQS" setextendeddbtype "${TMP_PATH}/aln_expanded" --extended-dbtype 2 ${VERBOSITY}
        fi
        INTERMEDIATE="${TMP_PATH}/aln_expanded"
    fi
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" "${ALIGNMENT_ALGO}" "${QUERYDB}" "${TARGET_ALIGNMENT}${INDEXEXT}" "$INTERMEDIATE" "${RESULTS}" ${TMP} \
          || fail "Alignment died"
  fi

  # create profiles
  if [ $STEP -eq 0 ]; then
    if notExists "$TMP_PATH/profile_${STEP}.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" result2profile "${QUERYDB}" "${TARGET_ALIGNMENT}${INDEXEXT}" "$TMP_PATH/aln_${STEP}" "$TMP_PATH/profile_${STEP}" ${VERBOSITY_THREADS_PAR} ${PROFILE_EVAL} \
                || fail "Create profile died"
    fi
  fi
	QUERYDB="$TMP_PATH/profile_${STEP}"
	STEP=$((STEP+1))
done


if [ -n "$REMOVE_TMP" ]; then
    STEP=0
    while [ "${STEP}" -lt "${NUM_IT}" ]; do
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${TMP_PATH}/pref_${STEP}" ${VERBOSITY}
      if [ $STEP -eq 0 ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/aln_${STEP}" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/profile_${STEP}" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/profile_${STEP}_ss" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/profile_${STEP}_h" ${VERBOSITY}
      fi
      STEP=$((STEP+1))
    done
    if [ -n "${EXPAND}" ]; then
       # shellcheck disable=SC2086
       "$MMSEQS" rmdb "${TMP_PATH}/aln_expanded" ${VERBOSITY}
    fi
    rm -f "$TMP_PATH/structureprofile.sh"
fi

