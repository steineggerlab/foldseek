#!/bin/sh
FLAGS="$(grep -m 1 '^flags' /proc/cpuinfo)"
case "${FLAGS}" in
  *avx2*)
    exec /usr/local/bin/foldseek_avx2 "$@"
    ;;
  *)
    exec /usr/local/bin/foldseek_sse41 "$@"
    ;;
esac
