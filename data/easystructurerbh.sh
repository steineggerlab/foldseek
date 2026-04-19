#!/bin/sh -e
# Foldseek-specific wrapper around easyrbh.sh that adds StrucTTY viewer support

# Run the original easyrbh workflow
sh "${TMP_PATH}/easyrbh.sh" "$@"

# View results with StrucTTY
if [ -n "${VIEW_RESULTS}" ] || [ -n "${STRUCTTY_PATH}" ]; then
    # Determine query structure file
    if [ -f "$1" ] && [ ! -f "${1}.dbtype" ]; then
        VIEWER_QUERY="$1"
    elif [ -n "${QUERY_INPUT}" ] && [ -f "${QUERY_INPUT}" ]; then
        VIEWER_QUERY="${QUERY_INPUT}"
    else
        VIEWER_QUERY=""
    fi

    # Determine target structure DB
    if [ -f "${TARGET}.dbtype" ]; then
        VIEWER_TARGET="${TARGET}"
    else
        VIEWER_TARGET=""
    fi

    sh "${TMP_PATH}/structty_viewer.sh" "${STRUCTTY_PATH}" "${VIEWER_QUERY}" "${RESULTS}" "${VIEWER_TARGET}"
fi
