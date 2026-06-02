#!/bin/sh -e
# Foldseek-specific wrapper around easyrbh.sh

# Cleanup-only re-invocation (after StrucTTY launch). Mirrors easyrbh.sh's tmp
# cleanup block so the tmp DBs survive until the viewer closes (D10). easyrbh.sh
# is an upstream (mmseqs) file and is intentionally not modified here.
if [ -n "${CLEANUP_ONLY}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result" ${VERBOSITY}
    if [ -z "${LEAVE_INPUT}" ]; then
        if [ -f "${TMP_PATH}/target" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_h" ${VERBOSITY}
        fi
        if [ -f "${TMP_PATH}/target_pad" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_pad" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_pad_h" ${VERBOSITY}
        fi
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_h" ${VERBOSITY}
        if [ -f "${TMP_PATH}/query_pad" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/query_pad" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/query_pad_h" ${VERBOSITY}
        fi
    fi
    rm -rf "${TMP_PATH}/rbh_tmp"
    rm -f "${TMP_PATH}/easyrbh.sh" "${TMP_PATH}/easystructurerbh.sh"
    exit 0
fi

# Run the original easyrbh workflow
sh "${TMP_PATH}/easyrbh.sh" "$@"
