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

# Merge two databases using soft-linking
# $1: input db1
# $2: input db2
# $3: output db
buildMergedDb() {
    # combine seq dbs
    MAXOFFSET=$(awk '($2+$3) > max{max=$2+$3}END{print max}' "${2}.index")
    awk -v OFFSET="${MAXOFFSET}" 'FNR==NR{print $0; next}{print $1"\t"$2+OFFSET"\t"$3}' "${2}.index" \
         "${1}.index" > "${3}.index"
    ln -s "$(abspath "${2}")" "${3}.0"
    ln -s "$(abspath "${1}")" "${3}.1"
    cp "${2}.dbtype" "${3}.dbtype"
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
      $RUNNER "$MMSEQS" structurerescorediagonal "${INPUT}" "${INPUT}" "${TMP_PATH}/pref" "${TMP_PATH}/pref_rescore1" ${STRUCTURERESCOREDIAGONAL_PAR} \
          || fail "Rescore with hamming distance step died"
  fi

  if notExists "${TMP_PATH}/pre_clust.dbtype"; then
      # shellcheck disable=SC2086,SC2153
      "$MMSEQS" clust "$INPUT" "${TMP_PATH}/pref_rescore1" "${TMP_PATH}/pre_clust" ${CLUSTER_PAR} \
          || fail "Pre-clustering step died"
  fi

  awk '{ print $1 }' "${TMP_PATH}/pre_clust.index" > "${TMP_PATH}/order_redundancy"

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
  if notExists "${TMP_PATH}/aln.linclust.dbtype"; then
      # shellcheck disable=SC2086
      $RUNNER "$MMSEQS" $ALIGNMENT_ALGO "${INPUT}${ALN_EXTENSION}" \
              "${INPUT}${ALN_EXTENSION}" "${TMP_PATH}/pref_filter2" \
              "${TMP_PATH}/aln.linclust" ${ALIGNMENT_PAR} || fail "Alignment step died"
  fi

  if notExists "${TMP_PATH}/pre_clustered_seqs.dbtype"; then
      # shellcheck disable=SC2086
      "$MMSEQS" createsubdb "${TMP_PATH}/order_redundancy" "${INPUT}" "${TMP_PATH}/pre_clustered_seqs" ${VERBOSITY} --subdb-mode 1 \
          || fail "Createsubdb pre_clustered_seqs step died"
  fi

  # 4. Clustering using greedy set cover.
  if notExists "${TMP_PATH}/clust.linclust.dbtype"; then
      # shellcheck disable=SC2086,SC2153
      "$MMSEQS" clust "${TMP_PATH}/pre_clustered_seqs" "${TMP_PATH}/aln.linclust" "${TMP_PATH}/clust.linclust" ${CLUSTER_PAR} \
          || fail "Clustering step died"
  fi

  if notExists "${TMP_PATH}/clu_redundancy.dbtype"; then
      # shellcheck disable=SC2086
      if [ "${RUN_ITERATIVE}" = "1" ]; then
         "$MMSEQS" mergeclusters "$SOURCE" "${TMP_PATH}/clu_redundancy" "${TMP_PATH}/pre_clust" "${TMP_PATH}/clust.linclust" $MERGECLU_PAR \
            || fail "mergeclusters died"
      else
         "$MMSEQS" mergeclusters "$SOURCE" "$2" "${TMP_PATH}/pre_clust" "${TMP_PATH}/clust.linclust" $MERGECLU_PAR \
            || fail "mergeclusters died"
      fi
  fi
fi

if [ "${RUN_ITERATIVE}" = "1" ]; then
  if [ "${RUN_LINCLUST}" = "1" ]; then
      if notExists "${TMP_PATH}/input_step_redundancy_ss.dbtype"; then
         # shellcheck disable=SC2086
         "$MMSEQS" createsubdb "${TMP_PATH}/clu_redundancy" "${INPUT}" "${TMP_PATH}/input_step_redundancy" ${VERBOSITY} --subdb-mode 1 \
            || fail "createsubdb died"
      fi
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
          $RUNNER "$MMSEQS" prefilter "${INPUT}_ss" "${INPUT}_ss" "${TMP_PATH}/pref_step$STEP" ${TMP} \
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
              "$MMSEQS" createsubdb "${TMP_PATH}/clu_step$STEP" "${INPUT}" "${NEXTINPUT}"  ${VERBOSITY} --subdb-mode 1 \
                                || fail "Order step $STEP died"
          fi
      fi

    INPUT="$NEXTINPUT"
    STEP=$((STEP+1))
  done
fi


if [ -n "$REASSIGN" ]; then
    STEP=$((STEP-1))
    PARAM=ALIGNMENT${STEP}_PAR
    eval ALIGNMENT_PAR="\$$PARAM"
    # align to cluster sequences
    if notExists "${TMP_PATH}/aln.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" $ALIGNMENT_ALGO "${SOURCE}${ALN_EXTENTION}" "${SOURCE}${ALN_EXTENTION}" "${TMP_PATH}/clu" "${TMP_PATH}/aln" ${ALIGNMENT_PAR} \
                || fail "Alignment step $STEP died"
    fi
    # create file of cluster that do not align based on given criteria
    if notExists "${TMP_PATH}/clu_not_accepted.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" subtractdbs "${TMP_PATH}/clu" "${TMP_PATH}/aln" "${TMP_PATH}/clu_not_accepted" --e-profile 100000000 -e 100000000 ${THREADSANDCOMPRESS} \
                 || fail "subtractdbs1 reassign died"
    fi
    if notExists "${TMP_PATH}/clu_not_accepted_swap.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" swapdb "${TMP_PATH}/clu_not_accepted" "${TMP_PATH}/clu_not_accepted_swap" ${THREADSANDCOMPRESS} \
                 || fail "swapdb1 reassign died"
    fi
    # short circuit if nothing can be reassigned
    if [ ! -s "${TMP_PATH}/clu_not_accepted_swap.index" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" mvdb "${TMP_PATH}/clu" "$2" ${VERBOSITY}
        if [ -n "$REMOVE_TMP" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/clu_not_accepted_swap" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/clu_not_accepted" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY}
        fi
    else
        # create file of cluster that do align based on given criteria
        if notExists "${TMP_PATH}/clu_accepted.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" subtractdbs "${TMP_PATH}/clu" "${TMP_PATH}/clu_not_accepted" "${TMP_PATH}/clu_accepted" --e-profile 100000000 -e 100000000 ${THREADSANDCOMPRESS} \
                     || fail "subtractdbs2 reassign died"
        fi
        # create sequences database that were wrong assigned
        if notExists "${TMP_PATH}/seq_wrong_assigned.dbtype"; then
            if notExists "${TMP_PATH}/seq_wrong_assigned_ss.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" createsubdb "${TMP_PATH}/clu_not_accepted_swap" "$SOURCE" "${TMP_PATH}/seq_wrong_assigned"  ${VERBOSITY} \
                     || fail "createsubdb1 reassign died"
            fi
        fi

        # build seed sequences
        if notExists "${TMP_PATH}/seq_seeds.dbtype"; then
            if notExists "${TMP_PATH}/seq_seeds_ss.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" createsubdb "${TMP_PATH}/clu" "$SOURCE" "${TMP_PATH}/seq_seeds" ${VERBOSITY}   \
                     || fail "createsubdb2 reassign died"
            fi

        fi
        PARAM=PREFILTER${STEP}_PAR
        eval PREFILTER_PAR="\$$PARAM"
        # try to find best matching centroid sequences for prev. wrong assigned sequences
        if notExists "${TMP_PATH}/seq_wrong_assigned_pref.dbtype"; then
            if notExists "${TMP_PATH}/seq_seeds.merged.dbtype"; then
                buildMergedDb "${TMP_PATH}/seq_wrong_assigned" "${TMP_PATH}/seq_seeds" "${TMP_PATH}/seq_seeds.merged"
                buildMergedDb "${TMP_PATH}/seq_wrong_assigned_ss" "${TMP_PATH}/seq_seeds_ss" "${TMP_PATH}/seq_seeds.merged_ss"
                if exists "${TMP_PATH}/seq_seeds_ca.dbtype"; then
                     buildMergedDb "${TMP_PATH}/seq_wrong_assigned_ca" "${TMP_PATH}/seq_seeds_ca" "${TMP_PATH}/seq_seeds.merged_ca"
                fi
            fi
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" prefilter "${TMP_PATH}/seq_wrong_assigned_ss" "${TMP_PATH}/seq_seeds.merged_ss" "${TMP_PATH}/seq_wrong_assigned_pref" ${PREFILTER_REASSIGN_PAR} \
                     || fail "Prefilter reassign died"
        fi
        if notExists "${TMP_PATH}/seq_wrong_assigned_pref_swaped.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" swapdb "${TMP_PATH}/seq_wrong_assigned_pref" "${TMP_PATH}/seq_wrong_assigned_pref_swaped" ${THREADSANDCOMPRESS} \
                     || fail "swapdb2 reassign died"
        fi
        if notExists "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln.dbtype"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" $ALIGNMENT_ALGO "${TMP_PATH}/seq_seeds.merged${ALN_EXTENTION}" "${TMP_PATH}/seq_wrong_assigned${ALN_EXTENTION}" \
                                              "${TMP_PATH}/seq_wrong_assigned_pref_swaped" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln" ${ALIGNMENT_REASSIGN_PAR} \
                     || fail "align2 reassign died"
        fi

        if notExists "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_ocol.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" filterdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_ocol" --trim-to-one-column ${THREADSANDCOMPRESS} \
                        || fail "filterdb2 reassign died"
        fi

        if notExists "${TMP_PATH}/clu_accepted_plus_wrong.dbtype"; then
            # combine clusters
            # shellcheck disable=SC2086
            "$MMSEQS" mergedbs "${TMP_PATH}/seq_seeds.merged" "${TMP_PATH}/clu_accepted_plus_wrong" "${TMP_PATH}/clu_accepted" \
                            "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_ocol" ${MERGEDBS_PAR} \
                                 || fail "mergedbs reassign died"
        fi

        if notExists "${TMP_PATH}/missing.single.seqs.db.dbtype"; then
             awk 'FNR==NR{if($3 > 1){ f[$1]=1; }next} !($1 in f){print $1"\t"$1}' "${TMP_PATH}/clu_accepted_plus_wrong.index" "${SOURCE}.index" > "${TMP_PATH}/missing.single.seqs"
            # shellcheck disable=SC2086
            "$MMSEQS" tsv2db "${TMP_PATH}/missing.single.seqs" "${TMP_PATH}/missing.single.seqs.db" --output-dbtype 6 ${VERBCOMPRESS} \
                                || fail "tsv2db reassign died"
        fi

        if notExists "${TMP_PATH}/clu_accepted_plus_wrong_plus_single.dbtype"; then
            # combine clusters
            # shellcheck disable=SC2086
            "$MMSEQS" mergedbs "${SOURCE}" "${TMP_PATH}/clu_accepted_plus_wrong_plus_single" "${TMP_PATH}/clu_accepted_plus_wrong" \
                            "${TMP_PATH}/missing.single.seqs.db" ${MERGEDBS_PAR} \
                                 || fail "mergedbs2 reassign died"
        fi

        PARAM=CLUSTER${STEP}_PAR
        eval TMP="\$$PARAM"
        # shellcheck disable=SC2086
        "$MMSEQS" clust "${SOURCE}" "${TMP_PATH}/clu_accepted_plus_wrong_plus_single" "${2}" ${TMP} \
                || fail "Clustering step $STEP died"

        if [ -n "$REMOVE_TMP" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/clu_not_accepted" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/clu_accepted" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/clu_not_accepted_swap" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/seq_seeds" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/seq_seeds.merged" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned_pref" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_ocol" ${VERBOSITY}
            rm -f "${TMP_PATH}/missing.single.seqs"
            rm -f "${TMP_PATH}/clu_accepted_plus_wrong.tsv"
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/missing.single.seqs.db" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/clu_accepted_plus_wrong" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/clu_accepted_plus_wrong_plus_single" ${VERBOSITY}
        fi
    fi
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
      "$MMSEQS" rmdb "${TMP_PATH}/aln.linclust" ${VERBOSITY}
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${TMP_PATH}/clust.linclust" ${VERBOSITY}
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${TMP_PATH}/pre_clustered_seqs" ${VERBOSITY}
    fi
    rm -f "${TMP_PATH}/clustering.sh"
fi


