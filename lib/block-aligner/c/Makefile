CC=clang
CFLAGS=-Ltarget/release -lblock_aligner -I.

.PHONY: all block_aligner example align_prefix

all: example align_prefix

block_aligner:
	cargo build --release --features simd_avx2 --offline
	cbindgen --config cbindgen.toml --crate block-aligner --output block_aligner.h --quiet .

example: block_aligner example.c
	$(CC) $(CFLAGS) -o example example.c

align_prefix: block_aligner align_prefix.c
	$(CC) $(CFLAGS) -o align_prefix align_prefix.c
