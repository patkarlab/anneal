#!/bin/bash
#
# docker -> apptainer shim
#
# anneal calls Parabricks through a hardcoded `docker run`. Docker is not
# usable from PBS batch jobs on this cluster (permission denied on the daemon
# socket), but apptainer is. Rather than patch and rebuild the locked binary,
# this shim sits earlier on PATH and rewrites the call.
#
# The invocation anneal emits (src/barcode/extract.rs) is fixed:
#
#   docker run --rm --gpus all \
#       --volume WORKDIR:/workdir --volume REFDIR:/refdir \
#       --workdir /workdir IMAGE pbrun fq2bam ...
#
# which becomes:
#
#   apptainer exec --nv -B WORKDIR:/workdir -B REFDIR:/refdir \
#       --pwd /workdir SIF pbrun fq2bam ...
#
# Any docker subcommand other than `run` (image inspect, load) exits 0 without
# doing anything, since apptainer needs no image preloading.
#
# Install:
#   mkdir -p ~/bin && cp docker_apptainer_shim.sh ~/bin/docker && chmod +x ~/bin/docker
#   export PATH="$HOME/bin:$PATH"      # must come before the real docker
#
# Override the image with PARABRICKS_SIF.

set -euo pipefail

SIF="${PARABRICKS_SIF:-$HOME/pipelines/parabricks_4.3.1.sif}"

if [ "${1:-}" != "run" ]; then
    # image inspect / load / anything else: nothing to do
    exit 0
fi
shift

if [ ! -f "$SIF" ]; then
    echo "docker-shim: SIF not found: $SIF" >&2
    exit 1
fi

NV=""
BINDS=()
PWDARG=()
CMD=()

while [ $# -gt 0 ]; do
    case "$1" in
        --rm|-i|-t|-it)
            shift
            ;;
        --gpus)
            NV="--nv"
            shift 2
            ;;
        --volume|-v)
            BINDS+=("-B" "$2")
            shift 2
            ;;
        --workdir|-w)
            PWDARG=("--pwd" "$2")
            shift 2
            ;;
        --env|-e)
            export "${2?}"
            shift 2
            ;;
        --*=*)
            shift
            ;;
        -*)
            echo "docker-shim: unhandled docker flag '$1', refusing to guess" >&2
            exit 1
            ;;
        *)
            # first non-flag is the image; everything after is the command
            shift
            CMD=("$@")
            break
            ;;
    esac
done

if [ ${#CMD[@]} -eq 0 ]; then
    echo "docker-shim: no command found after image name" >&2
    exit 1
fi

echo "docker-shim: apptainer exec ${NV} ${BINDS[*]} ${PWDARG[*]} $SIF ${CMD[*]}" >&2

exec apptainer exec ${NV} "${BINDS[@]}" "${PWDARG[@]}" "$SIF" "${CMD[@]}"
