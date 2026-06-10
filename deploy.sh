#!/bin/bash
# deploy.sh - Build Anneal 0.1.0
#
# Usage:
#   bash deploy.sh          # CPU-only build
#   bash deploy.sh gpu      # GPU build (uses pre-compiled PTX)
set -euo pipefail

BUILD_MODE="${1:-cpu}"

echo "=== Anneal 0.1.0 Build ==="
echo "Build mode: $BUILD_MODE"
echo ""

if ! command -v cargo &> /dev/null; then
    echo "ERROR: cargo not found. Activate conda environment first."
    exit 1
fi

if [ "$BUILD_MODE" = "gpu" ]; then
    PTX="src/cuda/consensus_kernel.ptx"

    # Check if real PTX exists (not placeholder)
    if [ -f "$PTX" ] && ! grep -q "Placeholder" "$PTX" 2>/dev/null; then
        echo "Using existing compiled PTX: $PTX"
    else
        # Allow an explicit prebuilt PTX via env var (optional)
        if [ -n "${ANNEAL_PTX:-}" ] && [ -f "${ANNEAL_PTX}" ] && ! grep -q "Placeholder" "${ANNEAL_PTX}" 2>/dev/null; then
            echo "Copying PTX from ANNEAL_PTX: ${ANNEAL_PTX}"
            cp "${ANNEAL_PTX}" "$PTX"
        fi

        # Still placeholder? Try nvcc
        if grep -q "Placeholder" "$PTX" 2>/dev/null; then
            if command -v nvcc &> /dev/null; then
                ARCH="${CUDA_ARCH:-sm_86}"
                echo "Compiling CUDA kernel (arch=$ARCH)..."
                nvcc -ptx -arch="$ARCH" src/cuda/consensus_kernel.cu -o "$PTX"
            else
                echo "ERROR: No compiled PTX found and nvcc not available."
                echo "Either:"
                echo "  1. Point ANNEAL_PTX at a prebuilt PTX, e.g.:"
                echo "     ANNEAL_PTX=/path/to/consensus_kernel.ptx bash deploy.sh gpu"
                echo "  2. Or compile on a node with nvcc + g++:"
                echo "     nvcc -ptx -arch=sm_86 src/cuda/consensus_kernel.cu -o src/cuda/consensus_kernel.ptx"
                exit 1
            fi
        fi
    fi

    echo ""
    echo "=== Building Anneal (GPU) ==="
    RUSTFLAGS="-C linker=gcc" cargo build --release --features gpu 2>&1
else
    echo ""
    echo "=== Building Anneal (CPU) ==="
    RUSTFLAGS="-C linker=gcc" cargo build --release 2>&1
fi

BINARY="target/release/anneal"
if [ -f "$BINARY" ]; then
    echo ""
    echo "=== Build successful ==="
    echo "Binary: $(pwd)/$BINARY"
    echo "Size: $(du -h "$BINARY" | cut -f1)"
else
    echo "ERROR: Binary not found"
    exit 1
fi
