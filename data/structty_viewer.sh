#!/bin/sh -e
# structty_viewer.sh — StrucTTY viewer common helper
# Usage: structty_viewer.sh [STRUCTTY_DIR] [QUERY_FILE] [RESULTS_M8] [TARGET_STRUCT_DB]
#
# Arguments:
#   STRUCTTY_DIR      — Directory containing StrucTTY binary (from --structty or --structty-dir option)
#                       Empty string means search PATH for "StrucTTY"
#   QUERY_FILE        — Query structure file path (PDB/CIF, empty string if none)
#   RESULTS_M8        — Foldseek .m8 result file path (required)
#   TARGET_STRUCT_DB  — Target structure DB path (empty string if none)
#                       StrucTTY uses --db to read _ca DB directly for local hit loading

STRUCTTY_DIR="$1"
QUERY_FILE="$2"
RESULTS_M8="$3"
TARGET_STRUCT_DB="$4"

# Determine binary: explicit directory first, then PATH search
if [ -n "${STRUCTTY_DIR}" ]; then
    STRUCTTY_BIN="${STRUCTTY_DIR}/StrucTTY"
    if [ ! -x "${STRUCTTY_BIN}" ]; then
        echo "Error: StrucTTY binary not found or not executable at ${STRUCTTY_BIN}"
        exit 1
    fi
elif command -v "StrucTTY" > /dev/null 2>&1; then
    STRUCTTY_BIN="StrucTTY"
else
    echo "Warning: StrucTTY not found in PATH. Install StrucTTY or use --structty <dir> to specify the directory containing StrucTTY."
    exit 0
fi

# Build command
STRUCTTY_CMD="${STRUCTTY_BIN}"

if [ -n "${QUERY_FILE}" ]; then
    STRUCTTY_CMD="${STRUCTTY_CMD} \"${QUERY_FILE}\""
fi

STRUCTTY_CMD="${STRUCTTY_CMD} --foldseek \"${RESULTS_M8}\""

if [ -n "${TARGET_STRUCT_DB}" ]; then
    STRUCTTY_CMD="${STRUCTTY_CMD} --db \"${TARGET_STRUCT_DB}\""
fi

eval "${STRUCTTY_CMD}"
