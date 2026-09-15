#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
container_args=()
for arg in "$@"; do
    case "$arg" in
        "$repo_root"/*)
            container_args+=("/work/${arg#"$repo_root/"}")
            ;;
        *)
            container_args+=("$arg")
            ;;
    esac
done

docker run --rm --platform linux/amd64 \
    -e RNAMOTIFSCANX_PATH=/work/external/rmsx/RNAMotifScanX_src \
    -v "$repo_root:/work" \
    ubuntu:22.04 \
    /work/external/rmsx/RNAMotifScanX_src/scan "${container_args[@]}"
