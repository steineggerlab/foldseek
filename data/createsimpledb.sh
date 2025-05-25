#!/bin/sh -e

if [ -e "${IN}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" createstructsimpledb "${IN}" "${OUT}"  ${VERBOSITY_PAR}  \
        || fail "createstructsimpledb died"
fi

if [ -e "${IN}.source" ]; then
    cp "${IN}.source" "${OUT}.source"
    awk '{print $0"\t"$1}' "${IN}.source" > "${OUT}.lookup"
    # shellcheck disable=SC2086
    "$MMSEQS" tsv2db "${IN}.source" "${OUT}_h" ${VERBOSITY_PAR} --output-dbtype 12 \
        || fail "tsv2b died"
fi

if [ -e "${OUT}.sh" ]; then
  rm -f -- "${OUT}.sh"
fi
