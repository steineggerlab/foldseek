#!/bin/sh -e

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

if [ -e "${IN}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" createStructinterfacedb "${IN}" "${OUT}" ${CREATESTRUCTINTERFACEDB_PAR} \
        || fail "createStructinterfacedb died"
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${OUT}.index" "${IN}_h" "${OUT}_h" \
        || fail "createsubdb died"
    
    awk 'BEGIN{
        i=0
        } FNR==NR{
            name[$1]=$2; next
            } FNR%2==1{
                print $1"\tINT"i"_"name[$1]"\t"i; getline; print $1"\tINT"i"_"name[$1]"\t"i; i++
                }' "$(abspath "${IN}.lookup")" "${OUT}.index" >  "${OUT}.lookup"
    awk 'NR%2 == 1{
        n = split($2, a, "_"); 
        hi=a[1]"_"
        for (i = 2; i < n-1; i++) {
            hi= hi a[i]"_";
        } 
        hi = hi a[n-1]
        print $3"\t"hi 
        }' "${OUT}.lookup" > "${OUT}.source"
fi

if [ -e "${OUT}.sh" ]; then
    rm -f -- "${OUT}.sh"
fi