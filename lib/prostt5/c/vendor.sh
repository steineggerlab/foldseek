#!/bin/sh -e
cargo-vendor-filterer vendor --features cuda,cudnn --platform aarch64-unknown-linux-gnu --platform x86_64-apple-darwin --platform x86_64-unknown-linux-gnu --platform aarch64-apple-darwin
find vendor \( -type d \( -name benches -o -name examples -o -name test -o -name tests \) -exec rm -rf -- "{}" + \) -o \( -type f -name "*.md" ! -path "vendor/safetensors/README.md" ! -path "vendor/bindgen_cuda/README.md" -delete \)
./checksum.py vendor
