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

name_args=()
if [[ -n "${RMSX_CONTAINER_NAME:-}" ]]; then
    name_args+=(--name "$RMSX_CONTAINER_NAME")
fi

docker run --rm --platform linux/amd64 "${name_args[@]}" \
    -e RNAMOTIFSCANX_PATH=/work/external/rmsx/RNAMotifScanX_src \
    -v "$repo_root:/work" \
    ubuntu:22.04 \
    /work/external/rmsx/RNAMotifScanX_src/scan "${container_args[@]}"
