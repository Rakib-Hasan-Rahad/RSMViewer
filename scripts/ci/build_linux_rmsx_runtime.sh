#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
runtime_dir="$repo_root/rsmviewer/tools/rmsx_runtime"
src_dir="$runtime_dir/src/RNAMotifScanX_src"
linux_bin_dir="$runtime_dir/bin/linux-x86_64"
linux_build_dir="$runtime_dir/build/linux-x86_64-ci"
report_dir="$repo_root/dist/reports"
artifact_dir="$repo_root/dist/rmsx_runtime_linux_x86_64"
artifact_lib_dir="$artifact_dir/lib"

mkdir -p "$linux_bin_dir" "$linux_build_dir" "$report_dir" "$artifact_dir" "$artifact_lib_dir"

copy_mccore_if_present() {
  local mccore_root="$runtime_dir/src/external/mccore-install-linux-x86_64"
  local lib_path=""
  for candidate in \
    "$mccore_root/lib64/libmccore.so" \
    "$mccore_root/lib/libmccore.so"; do
    if [[ -f "$candidate" ]]; then
      lib_path="$candidate"
      break
    fi
  done
  if [[ -n "$lib_path" ]]; then
    cp "$lib_path" "$linux_bin_dir/libmccore.so"
    cp "$lib_path" "$artifact_lib_dir/libmccore.so"
  fi
}

patch_rpath_if_possible() {
  local bin_path="$1"
  if command -v patchelf >/dev/null 2>&1; then
    patchelf --set-rpath '$ORIGIN:$ORIGIN/lib' "$bin_path" || true
  fi
}

verify_binary_links() {
  local bin_path="$1"
  local report_path="$2"
  ldd "$bin_path" | tee "$report_path"
  if grep -q "not found" "$report_path"; then
    echo "[linux-runtime] ERROR: unresolved shared libraries in $bin_path" >&2
    return 1
  fi
}

bundle_runtime_libs() {
  local bin_path="$1"
  local deps
  deps="$(ldd "$bin_path" | awk '/=>/ {print $3} /^[[:space:]]*\// {print $1}' | sed '/^$/d' | sort -u)"

  while IFS= read -r dep; do
    [[ -z "$dep" ]] && continue
    local base
    base="$(basename "$dep")"
    case "$base" in
      libboost_*|libz.so*|libstdc++.so*|libgcc_s.so*|libmccore.so*)
        cp -n "$dep" "$artifact_lib_dir/$base" || true
        ;;
      *)
        ;;
    esac
  done <<< "$deps"
}

echo "[linux-runtime] Building required runtime pieces via setup_runtime.py"
python3 "$runtime_dir/setup_runtime.py" \
  --runtime-dir "$runtime_dir" \
  --build \
  --fetch-mccore \
  --fetch-mc-annotate \
  --skip-rnaview \
  --json > "$report_dir/setup_build.json"

echo "[linux-runtime] Building RNAVIEW for Linux"
python3 "$runtime_dir/setup_runtime.py" \
  --runtime-dir "$runtime_dir" \
  --rebuild-rnaview \
  --json > "$report_dir/setup_with_rnaview.json"

echo "[linux-runtime] Building cut and align with CMake"
cmake -S "$src_dir" -B "$linux_build_dir" -DCMAKE_BUILD_TYPE=Release -DBoost_NO_BOOST_CMAKE=ON
cmake --build "$linux_build_dir" --config Release --target cut align -- -j"$(nproc)"

copy_from_build() {
  local target_name="$1"
  local source_path=""
  if [[ -f "$linux_build_dir/$target_name" ]]; then
    source_path="$linux_build_dir/$target_name"
  elif [[ -f "$linux_build_dir/Release/$target_name" ]]; then
    source_path="$linux_build_dir/Release/$target_name"
  else
    echo "[linux-runtime] ERROR: built target not found: $target_name" >&2
    return 1
  fi
  cp "$source_path" "$linux_bin_dir/$target_name"
  chmod +x "$linux_bin_dir/$target_name"
}

copy_from_build cut
copy_from_build align

copy_mccore_if_present

patch_rpath_if_possible "$linux_bin_dir/scan"
patch_rpath_if_possible "$linux_bin_dir/cut"
patch_rpath_if_possible "$linux_bin_dir/align"
patch_rpath_if_possible "$linux_bin_dir/rnaview"
patch_rpath_if_possible "$linux_bin_dir/MC-Annotate"

echo "[linux-runtime] Validating runtime readiness"
python3 "$runtime_dir/setup_runtime.py" \
  --runtime-dir "$runtime_dir" \
  --json > "$report_dir/final_runtime_report.json"

python3 - "$report_dir/final_runtime_report.json" <<'PY'
import json
import sys
from pathlib import Path

report_path = Path(sys.argv[1])
report = json.loads(report_path.read_text(encoding="utf-8"))
if not report.get("ok", False):
    missing = ", ".join(report.get("missing", []))
    raise SystemExit(f"Runtime validation failed; missing: {missing}")
PY

required_bins=(scan cut align rnaview MC-Annotate)
rm -f "$report_dir/ldd_summary.txt"
for bin_name in "${required_bins[@]}"; do
  path="$linux_bin_dir/$bin_name"
  if [[ ! -x "$path" ]]; then
    echo "[linux-runtime] ERROR: missing executable $path" >&2
    exit 1
  fi
  file "$path"
  verify_binary_links "$path" "$report_dir/ldd_${bin_name}.txt"
  cat "$report_dir/ldd_${bin_name}.txt" >> "$report_dir/ldd_summary.txt"
  printf '\n' >> "$report_dir/ldd_summary.txt"
  bundle_runtime_libs "$path"
done

cp "$linux_bin_dir/scan" "$artifact_dir/"
cp "$linux_bin_dir/cut" "$artifact_dir/"
cp "$linux_bin_dir/align" "$artifact_dir/"
cp "$linux_bin_dir/rnaview" "$artifact_dir/"
cp "$linux_bin_dir/MC-Annotate" "$artifact_dir/"
cp "$report_dir/final_runtime_report.json" "$artifact_dir/runtime_report.json"
cp "$report_dir/ldd_summary.txt" "$artifact_dir/ldd_summary.txt"

if [[ -f "$linux_bin_dir/libmccore.so" ]]; then
  cp "$linux_bin_dir/libmccore.so" "$artifact_dir/"
fi

(
  cd "$repo_root/dist"
  tar -czf rmsx_runtime_linux_x86_64.tar.gz rmsx_runtime_linux_x86_64
)

echo "[linux-runtime] Done: dist/rmsx_runtime_linux_x86_64.tar.gz"
