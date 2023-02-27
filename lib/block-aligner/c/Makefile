CC=clang
CFLAGS=-L../target/release -lblock_aligner -I.

.PHONY: all block_aligner example align_prefix

all: example align_prefix

block_aligner:
	cargo build --release --features c_bindings,simd_avx2

example: block_aligner example.c
	$(CC) $(CFLAGS) -o example example.c

align_prefix: block_aligner align_prefix.c
	$(CC) $(CFLAGS) -o align_prefix align_prefix.c
