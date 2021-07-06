#!/bin/sh -e
# Assembler workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}
