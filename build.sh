#!/bin/bash
# build.sh -- build the anneal consensus binary (CPU build, the validated one)
#
#   bash build.sh            # builds target_cpu/release/anneal and runs the unit tests
#   bash build.sh --no-test  # build only
#
# The GPU cargo feature exists but is not used: consensus runs on the CPU in
# every validated configuration and stage 1 passes --no-gpu. config.sh points
# ANNEAL at target_cpu/release/anneal.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

command -v cargo >/dev/null 2>&1 || { echo "ERROR: cargo not found (source ~/.cargo/env)"; exit 1; }

export CARGO_TARGET_DIR=target_cpu
export RUSTFLAGS="-C linker=gcc"

echo "=== building anneal $(grep -m1 '^version' Cargo.toml | cut -d'"' -f2) (CPU) ==="
cargo build --release
if [ "${1:-}" != "--no-test" ]; then
    echo "=== unit tests ==="
    cargo test --release 2>&1 | tail -3
fi
echo "=== binary ==="
ls -l target_cpu/release/anneal
target_cpu/release/anneal --version 2>/dev/null || true
