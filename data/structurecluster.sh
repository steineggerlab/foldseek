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


# check number of input variables
[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[   -f "$2.dbtype" ] && echo "$2.dbtype exists already!" && exit 1;
[ ! -d "$3" ] && echo "tmp directory $3 not found!" && mkdir -p "$3";

INPUT="$1"
TMP_PATH="$3"
SOURCE="$INPUT"

if [ "${RUN_LINCLUST}" = "1" ]; then

  # 1. Finding exact $k$-mer matches.
  if notExists "${TMP_PATH}/pref.dbtype"; then
      # shellcheck disable=SC2086
      $RUNNER "$MMSEQS" kmermatcher "${INPUT}_ss" "${TMP_PATH}/pref" ${KMERMATCHER_PAR} \
          || fail "kmermatcher died"
  fi

  # 2. Hamming distance pre-clustering
  if notExists "${TMP_PATH}/pref_rescore1.dbtype"; then
      # shellcheck disable=SC2086
      $RUNNER "$MMSEQS" rescorediagonal "${INPUT}_ss" "${INPUT}_ss" "${TMP_PATH}/pref" "${TMP_PATH}/pref_rescore1" ${HAMMING_PAR} \
          || fail "Rescore with hamming distance step died"
  fi

  if notExists "${TMP_PATH}/pre_clust.dbtype"; then
      # shellcheck disable=SC2086,SC2153
      "$MMSEQS" clust "$INPUT" "${TMP_PATH}/pref_rescore1" "${TMP_PATH}/pre_clust" ${CLUSTER_PAR} \
          || fail "Pre-clustering step died"
  fi

  awk '{ print $1 }' "${TMP_PATH}/pre_clust.index" > "${TMP_PATH}/order_redundancy"
  if notExists "${TMP_PATH}/input_step_redundancy.dbtype"; then
      # shellcheck disable=SC2086
      "$MMSEQS" createsubdb "${TMP_PATH}/order_redundancy" "${INPUT}_ss" "${TMP_PATH}/input_step_redundancy_ss" ${VERBOSITY} --subdb-mode 1 \
          || fail "createsubdb step died"
      # shellcheck disable=SC2086
      "$MMSEQS" createsubdb "${TMP_PATH}/order_redundancy" "${INPUT}_ca" "${TMP_PATH}/input_step_redundancy_ca" ${VERBOSITY} --subdb-mode 1 \
          || fail "createsubdb step died"
      # shellcheck disable=SC2086
      "$MMSEQS" createsubdb "${TMP_PATH}/order_redundancy" "${INPUT}" "${TMP_PATH}/input_step_redundancy" ${VERBOSITY} --subdb-mode 1 \
          || fail "createsubdb step died"
  fi

  if notExists "${TMP_PATH}/pref_filter1.dbtype"; then
      # shellcheck disable=SC2086
      "$MMSEQS" createsubdb "${TMP_PATH}/order_redundancy" "${TMP_PATH}/pref" "${TMP_PATH}/pref_filter1" ${VERBOSITY} --subdb-mode 1 \
          || fail "Createsubdb step died"
  fi

  if notExists "${TMP_PATH}/pref_filter2.dbtype"; then
      # shellcheck disable=SC2086
      "$MMSEQS" filterdb "${TMP_PATH}/pref_filter1" "${TMP_PATH}/pref_filter2" --filter-file "${TMP_PATH}/order_redundancy" ${VERBOSITYANDCOMPRESS} \
          || fail "Filterdb step died"
  fi

  # 3. Local gapped sequence alignment.
  if notExists "${TMP_PATH}/aln.dbtype"; then
      # shellcheck disable=SC2086
      $RUNNER "$MMSEQS" $ALIGNMENT_ALGO "${TMP_PATH}/input_step_redundancy${ALN_EXTENSION}" \
              "${TMP_PATH}/input_step_redundancy${ALN_EXTENSION}" "${TMP_PATH}/pref_filter2" \
              "${TMP_PATH}/aln" ${ALIGNMENT_PAR} || fail "Alignment step died"
  fi
  # 4. Clustering using greedy set cover.
  if notExists "${TMP_PATH}/clust.dbtype"; then
      # shellcheck disable=SC2086,SC2153
      "$MMSEQS" clust "${TMP_PATH}/input_step_redundancy_ss" "${TMP_PATH}/aln" "${TMP_PATH}/clust" ${CLUSTER_PAR} \
          || fail "Clustering step died"
  fi
  if notExists "${TMP_PATH}/clu.dbtype"; then
      # shellcheck disable=SC2086
      if [ "${RUN_ITERATIVE}" = "1" ]; then
         "$MMSEQS" mergeclusters "$SOURCE" "${TMP_PATH}/clu_redundancy" "${TMP_PATH}/pre_clust" "${TMP_PATH}/clust" $MERGECLU_PAR \
            || fail "mergeclusters died"
      else
         "$MMSEQS" mergeclusters "$SOURCE" "$2" "${TMP_PATH}/pre_clust" "${TMP_PATH}/clust" $MERGECLU_PAR \
            || fail "mergeclusters died"
      fi
  fi
fi

if [ "${RUN_ITERATIVE}" = "1" ]; then
  if [ "${RUN_LINCLUST}" = "1" ]; then
      # shellcheck disable=SC2086
      "$MMSEQS" createsubdb "${TMP_PATH}/clu_redundancy" "${INPUT}_ss" "${TMP_PATH}/input_step_redundancy_ss" ${VERBOSITY} --subdb-mode 1 \
          || fail "createsubdb died"
      # shellcheck disable=SC2086
      "$MMSEQS" createsubdb "${TMP_PATH}/clu_redundancy" "${INPUT}_ca" "${TMP_PATH}/input_step_redundancy_ca" ${VERBOSITY} --subdb-mode 1 \
                || fail "createsubdb died"
      # shellcheck disable=SC2086
      "$MMSEQS" createsubdb "${TMP_PATH}/clu_redundancy" "${INPUT}" "${TMP_PATH}/input_step_redundancy" ${VERBOSITY} --subdb-mode 1 \
              || fail "createsubdb died"
      INPUT="${TMP_PATH}/input_step_redundancy"
  fi
  STEP=0
  STEPS=${STEPS:-1}
  CLUSTER_STR=""
  while [ "$STEP" -lt "$STEPS" ]; do
      PARAM=PREFILTER${STEP}_PAR
      eval TMP="\$$PARAM"
      if notExists "${TMP_PATH}/pref_step$STEP.dbtype"; then
           # shellcheck disable=SC2086
          $RUNNER "$MMSEQS" prefilter "${INPUT}" "${INPUT}" "${TMP_PATH}/pref_step$STEP" ${TMP} \
              || fail "Prefilter step $STEP died"
      fi
      PARAM=ALIGNMENT${STEP}_PAR
      eval TMP="\$$PARAM"
      if notExists "${TMP_PATH}/aln_step$STEP.dbtype"; then
          # shellcheck disable=SC2086
          $RUNNER "$MMSEQS" $ALIGNMENT_ALGO "${INPUT}${ALN_EXTENTION}" "${INPUT}${ALN_EXTENTION}" "${TMP_PATH}/pref_step$STEP" "${TMP_PATH}/aln_step$STEP" ${ALIGNMENT_PAR} \
              || fail "Alignment step $STEP died"
      fi
      PARAM=CLUSTER${STEP}_PAR
      eval TMP="\$$PARAM"
      if notExists "${TMP_PATH}/clu_step$STEP.dbtype"; then
           # shellcheck disable=SC2086
          "$MMSEQS" clust "${INPUT}" "${TMP_PATH}/aln_step$STEP" "${TMP_PATH}/clu_step$STEP" ${TMP} \
              || fail "Clustering step $STEP died"
      fi

      # FIXME: This won't work if paths contain spaces
      CLUSTER_STR="${CLUSTER_STR} ${TMP_PATH}/clu_step$STEP"
      NEXTINPUT="${TMP_PATH}/input_step$((STEP+1))"
      if [ "$STEP" -eq "$((STEPS-1))" ]; then
         if [ -n "$REASSIGN" ]; then
            if notExists "${TMP_PATH}/clu.dbtype"; then
              # shellcheck disable=SC2086
              "$MMSEQS" mergeclusters "$SOURCE" "${TMP_PATH}/clu" "${TMP_PATH}/clu_redundancy" ${CLUSTER_STR} \
              || fail "Merging of clusters has died"
            fi
         else
              # shellcheck disable=SC2086
              "$MMSEQS" mergeclusters "$SOURCE" "$2" "${TMP_PATH}/clu_redundancy" ${CLUSTER_STR} $MERGECLU_PAR \
              || fail "Merging of clusters has died"
         fi
      else
          if notExists "$NEXTINPUT.dbtype"; then
              # shellcheck disable=SC2086
              "$MMSEQS" createsubdb "${TMP_PATH}/clu_step$STEP" "${INPUT}_ss" "${NEXTINPUT}_ss" ${VERBOSITY} --subdb-mode 1 \
                  || fail "Order step $STEP died"
              # shellcheck disable=SC2086
              "$MMSEQS" createsubdb "${TMP_PATH}/clu_step$STEP" "${INPUT}_ca" "${NEXTINPUT}_ca" ${VERBOSITY} --subdb-mode 1 \
                                || fail "Order step $STEP died"
              # shellcheck disable=SC2086
              "$MMSEQS" createsubdb "${TMP_PATH}/clu_step$STEP" "${INPUT}" "${NEXTINPUT}" ${VERBOSITY} --subdb-mode 1 \
                                || fail "Order step $STEP died"
          fi
      fi

    INPUT="$NEXTINPUT"
    STEP=$((STEP+1))
  done
fi

if [ -n "$REMOVE_TMP" ]; then
    if [ "${RUN_ITERATIVE}" = "1" ]; then
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${TMP_PATH}/clu_redundancy" ${VERBOSITY}
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${TMP_PATH}/input_step_redundancy" ${VERBOSITY}
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${TMP_PATH}/input_step_redundancy_h" ${VERBOSITY}
      STEP=0
      while [ "$STEP" -lt "$STEPS" ]; do
          # shellcheck disable=SC2086
          "$MMSEQS" rmdb "${TMP_PATH}/pref_step$STEP" ${VERBOSITY}
          # shellcheck disable=SC2086
          "$MMSEQS" rmdb "${TMP_PATH}/aln_step$STEP" ${VERBOSITY}
          # shellcheck disable=SC2086
          "$MMSEQS" rmdb "${TMP_PATH}/clu_step$STEP" ${VERBOSITY}
          STEP=$((STEP+1))
      done

      STEP=1
      while [ "$STEP" -lt "$STEPS" ]; do
          # shellcheck disable=SC2086
          "$MMSEQS" rmdb "${TMP_PATH}/input_step$STEP" ${VERBOSITY}
          # shellcheck disable=SC2086
          "$MMSEQS" rmdb "${TMP_PATH}/input_step${STEP}_h" ${VERBOSITY}
          STEP=$((STEP+1))
      done
    fi
    if [ "${RUN_LINCLUST}" = "1" ]; then
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${TMP_PATH}/pref_filter1" ${VERBOSITY}
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${TMP_PATH}/pref" ${VERBOSITY}
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${TMP_PATH}/pref_rescore1" ${VERBOSITY}
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${TMP_PATH}/pre_clust" ${VERBOSITY}
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${TMP_PATH}/input_step_redundancy" ${VERBOSITY}
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${TMP_PATH}/input_step_redundancy_h" ${VERBOSITY}
      rm -f "${TMP_PATH}/order_redundancy"
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${TMP_PATH}/pref_filter2" ${VERBOSITY}
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY}
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${TMP_PATH}/clust" ${VERBOSITY}
    fi
    rm -f "${TMP_PATH}/clustering.sh"
fi

