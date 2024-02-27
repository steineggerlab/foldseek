#!/bin/sh -e
#TODO: maybe change file name into filtercomplex.sh
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

# Shift initial DB to complexDB using soft-linking
# $1: input db
# $2: output db
buildCmplDb() {
    touch "${2}"
    awk -F"\t" 'BEGIN {OFFSET=0}
        FNR==NR{chain_len[$1]=$3;next}
        {
            if (!($3 in off_arr)) {
                off_arr[$3]=OFFSET
            }
            cmpl_len[$3]+=chain_len[$1];OFFSET+=chain_len[$1]
        }
        END {
            for (cmpl in off_arr) {
                print cmpl"\t"off_arr[cmpl]"\t"cmpl_len[cmpl]
            }
        }' "${1}.index" "${1}.lookup" > "${2}.index"
    ln -s "$(abspath "${1}")" "${2}.0"
    cp "${1}.dbtype" "${2}.dbtype"
}

# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <alignmentDB> <clustDB>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[ ! -f "$3.dbtype" ] && echo "$3.dbtype not found!" && exit 1;
[   -f "$4.dbtype" ] && echo "$4.dbtype exists already!" && exit 1;

if notExists "$4"; then
    # shellcheck disable=SC2086
    $MMSEQS filtercomplex "$1" "$2" "$3" "$4" ${FILTERCOMPLEX_PAR} \
        || fail "FilterComplex died"
fi

if notExists "${CMPLDB_PATH}/cmpl_db.dbtype"; then
    buildCmplDb "${SOURCE}" "${CMPLDB_PATH}/cmpl_db"
fi

# DONE : remove tmp -> No TMP file generated