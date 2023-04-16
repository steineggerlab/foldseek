#!/bin/sh -e
# Assembler workflow script
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

if notExists "${RESULTDB}_seq.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" cpdb "${INPUTDB}" "${RESULTDB}_seq" ${VERBOSITY}
fi

if notExists "${RESULTDB}_seq_h.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" cpdb "${INPUTDB}_h" "${RESULTDB}_seq_h" ${VERBOSITY}
fi

if notExists "${RESULTDB}_seq_ss.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" cpdb "${INPUTDB}_ss" "${RESULTDB}_seq_ss" ${VERBOSITY}
fi

if notExists "${RESULTDB}_seq_ca.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" cpdb "${INPUTDB}_ca" "${RESULTDB}_seq_ca" ${VERBOSITY}
fi

if notExists "${RESULTDB}_clu.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" cpdb "${CLUSTERDB}" "${RESULTDB}_clu" ${VERBOSITY}
fi

if notExists "${RESULTDB}_ca.dbtype"; then
  # shellcheck disable=SC2086
  "$MMSEQS" createsubdb "${CLUSTERDB}" "${RESULTDB}_seq_ca" "${RESULTDB}_ca" ${VERBOSITY}
fi

if notExists "${RESULTDB}_aln.dbtype"; then
  # shellcheck disable=SC2086
  "$MMSEQS" structurealign  "${RESULTDB}_seq" "${RESULTDB}_seq" "${CLUSTERDB}" "${RESULTDB}_aln" -a -e 0.1 --sort-by-structure-bits 0 ${VERBOSITYANDTHREADS}
fi

if notExists "${RESULTDB}_profile.dbtype"; then
  # shellcheck disable=SC2086
  "$MMSEQS" result2profile  "${RESULTDB}_seq" "${RESULTDB}_seq" "${RESULTDB}_aln" "${RESULTDB}_profile" ${PROFILE_PAR}
fi

if notExists "${RESULTDB}.dbtype"; then
  # shellcheck disable=SC2086
  "$MMSEQS" profile2consensus "${RESULTDB}_profile" "${RESULTDB}" ${VERBOSITYANDTHREADS}
  if [ -n "$REMOVE_TMP" ]; then
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${RESULTDB}_profile" ${VERBOSITY}
  fi
fi

if notExists "${RESULTDB}_profile_ss.dbtype"; then
  # shellcheck disable=SC2086
  "$MMSEQS" result2profile  "${RESULTDB}_seq_ss" "${RESULTDB}_seq_ss" "${RESULTDB}_aln" "${RESULTDB}_profile_ss" ${PROFILE_SS_PAR}
fi

if notExists "${RESULTDB}_ss.dbtype"; then
  # shellcheck disable=SC2086
  "$MMSEQS" profile2consensus "${RESULTDB}_profile_ss" "${RESULTDB}_ss" ${VERBOSITYANDTHREADS}
  if [ -n "$REMOVE_TMP" ]; then
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${RESULTDB}_profile_ss" ${VERBOSITY}
      # shellcheck disable=SC2086
      "$MMSEQS" rmdb "${RESULTDB}_aln" ${VERBOSITY}
  fi
fi

