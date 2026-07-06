#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

exists() {
	[ -f "$1" ]
}

notExists() {
	[ ! -f "$1" ]
}

# Remove the tmp DBs produced by this workflow. Defined as a function so it can be
# invoked both at the end of a normal run (REMOVE_TMP) and on a CLEANUP_ONLY re-run
# after the StrucTTY viewer closes (D10: launch must happen before cleanup).
do_cleanup() {
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/multimer_result" ${VERBOSITY}
    if [ -z "${LEAVE_INPUT}" ]; then
        if [ -f "${TMP_PATH}/target" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_h" ${VERBOSITY}
            if exists "${TMP_PATH}/target_ca.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/target_ca" ${VERBOSITY}
            fi
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_ss" ${VERBOSITY}
        fi
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_h" ${VERBOSITY}
        if exists "${TMP_PATH}/target_ca.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/query_ca" ${VERBOSITY}
        fi
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_ss" ${VERBOSITY}
    fi
    rm -rf "${TMP_PATH}/multimersearch_tmp"
    rm -f "${TMP_PATH}/easymultimersearch.sh"
}

# Cleanup-only re-invocation (after StrucTTY launch). Run cleanup and exit.
if [ -n "${CLEANUP_ONLY}" ]; then
    do_cleanup
    exit 0
fi

if notExists "${QUERY}.dbtype"; then
    if notExists "${TMP_PATH}/query"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createdb "${QUERY}" "${TMP_PATH}/query" ${CREATEDB_PAR} \
            || fail "query createdb died"
    fi
    QUERY="${TMP_PATH}/query"
fi

if notExists "${TARGET}.dbtype"; then
    if notExists "${TMP_PATH}/target"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createdb "${TARGET}" "${TMP_PATH}/target" ${CREATEDB_PAR} \
            || fail "target createdb died"
    fi
    TARGET="${TMP_PATH}/target"

    if [ -n "${GPU}" ]; then
        if notExists "${TMP_PATH}/target_pad"; then
            # shellcheck disable=SC2086
            "$MMSEQS" makepaddedseqdb "${TMP_PATH}/target" "${TMP_PATH}/target_pad" ${MAKEPADDEDSEQDB_PAR} \
                || fail "makepaddedseqdb died"
        fi
        TARGET="${TMP_PATH}/target_pad"
    fi
fi

if notExists "${TMP_PATH}/multimer_result.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" multimersearch "${QUERY}" "${TARGET}" "${TMP_PATH}/multimer_result" "${TMP_PATH}/multimersearch_tmp" ${MULTIMERSEARCH_PAR} \
    || fail "multimersearch died"
fi

# shellcheck disable=SC2086
"$MMSEQS" convertalis "${QUERY}" "${TARGET}" "${TMP_PATH}/multimer_result" "${OUTPUT}" ${CONVERT_PAR} \
    || fail "Convert Alignments died"

if [ -z "${NO_REPORT}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" createmultimerreport "${QUERY}" "${TARGET}" "${TMP_PATH}/multimer_result" "${OUTPUT}_report" ${REPORT_PAR} \
        || fail "createmultimerreport died"
fi

# View results with StrucTTY
if [ -n "${VIEW_RESULTS}" ]; then
    # Determine StrucTTY binary: prefer STRUCTTY_PATH env var, then PATH lookup
    if [ -n "${STRUCTTY_PATH}" ]; then
        if [ -x "${STRUCTTY_PATH}" ]; then
            STRUCTTY_BIN="${STRUCTTY_PATH}"
        else
            echo "Error: StrucTTY binary not found or not executable at ${STRUCTTY_PATH}"
            STRUCTTY_BIN=""
        fi
    elif command -v "StrucTTY" > /dev/null 2>&1; then
        STRUCTTY_BIN="StrucTTY"
    else
        STRUCTTY_BIN=""
    fi

    if [ -n "${STRUCTTY_BIN}" ]; then
        if notExists "${OUTPUT}_report"; then
            echo "Warning: --view-structty requires the multimer report (${OUTPUT}_report); skipping StrucTTY launch."
            echo "Results have been saved to: ${OUTPUT}"
        else
            UT_FILE="${TMP_PATH}/structty_ut.tsv"
            COMPLEX_M8="${TMP_PATH}/structty_complex.m8"

            # Extract per-complex U/T matrices (columns 7 and 8 of createmultimerreport output)
            awk -F'\t' 'BEGIN{OFS="\t"} { print NR, NR, $7, $8 }' \
                "${OUTPUT}_report" > "${UT_FILE}" \
                || fail "failed to extract u,t matrices for StrucTTY"

            # Generate complex-level m8: one row per complex pair (not per chain)
            # query=qComplexName($1), target=tComplexName($2), fident=qTMScore($5), rest=0
            awk -F'\t' 'BEGIN{OFS="\t"} { print $1, $2, $5, 0, 0, 0, 0, 0, 0, 0, 0, 0 }' \
                "${OUTPUT}_report" > "${COMPLEX_M8}" \
                || fail "failed to generate complex-level m8 for StrucTTY"

            STRUCTTY_CMD="${STRUCTTY_BIN}"
            if [ -n "${QUERY_INPUT}" ]; then
                STRUCTTY_CMD="${STRUCTTY_CMD} \"${QUERY_INPUT}\""
            fi
            if [ -n "${TARGET_INPUT}" ]; then
                STRUCTTY_CMD="${STRUCTTY_CMD} \"${TARGET_INPUT}\""
            fi
            STRUCTTY_CMD="${STRUCTTY_CMD} -ut \"${UT_FILE}\" --foldseek \"${COMPLEX_M8}\""
            if exists "${TARGET}.dbtype"; then
                STRUCTTY_CMD="${STRUCTTY_CMD} --db \"${TARGET}\""
            fi

            eval "${STRUCTTY_CMD}"
            rm -f "${UT_FILE}" "${COMPLEX_M8}"
        fi
    else
        echo "Warning: StrucTTY not found in PATH. Install StrucTTY or use --structty <path> to specify the binary location."
        echo "Results have been saved to: ${OUTPUT}"
    fi
fi

if [ -n "${REMOVE_TMP}" ]; then
    do_cleanup
fi