# RMSX Runtime Assets (Source 7)

This folder contains integrated runtime assets for Source 7 in RSMViewer.

## Bundled assets

- `queries/`: RNAMotifScanX motif query files (`*.struct`)
- `mat/`: scoring matrix files used by the RMSX toolchain
- `bin/<platform>/`: platform-specific executables (to be provided)

## Expected executables by platform

RSMViewer auto-detects these filenames in `bin/<platform>/`:

- RMSX executable: `RNAMotifScanX`, `scan`, `RNAMotifScanX.exe`, or `scan.exe`
- MC-Annotate executable: `MC-Annotate`, `mc-annotate`, `MC-Annotate.exe`, or `mc-annotate.exe`

Platform folders:

- `bin/macos-arm64`
- `bin/macos-x86_64`
- `bin/linux-x86_64`
- `bin/windows-x86_64`

## Source 7 runtime behavior

Source 7 runs in strict mode:

- no external config required
- no cached-result fallback
- fresh pipeline execution on each run
- no external CIF download during Source-7 runs

If binaries are missing for the current platform, Source 7 run fails fast with an explicit error.

## Build System

RNAMotifScanX C++ sources now include a cross-platform CMake project at:

- `src/RNAMotifScanX_src/CMakeLists.txt`

`rmv_rmsx setup` uses this CMake build flow to build `scan`/`scan.exe` on:

- Windows (MSVC)
- macOS (Apple Clang)
- Linux (GCC/Clang)

Required build dependencies:

- CMake
- C++ compiler toolchain for your platform
- Boost (program_options, thread, system, filesystem, iostreams)
- zlib

## Current Upstream Limitation

The source bundle currently ships MC-Annotate only as prebuilt Linux ELF binaries:

- `src/RNAMotifScanX_src/StructureAnnotation/MC-Annotate`
- `src/RNAMotifScanX_src/ThirdParty/MC-Annotate`

There is no MC-Annotate source code in this bundle to rebuild for Windows/macOS. As a result:

- Linux can use the shipped MC-Annotate binary.
- Windows/macOS require a native MC-Annotate executable to be supplied separately.

`rmv_rmsx_doctor` and `rmv_rmsx setup` report this explicitly when detected.

## Optional MC-Annotate Source Build

If you allow external source fetch during setup, you can ask setup script to clone/build MC-Annotate from:

- `https://github.com/major-lab/MC-Annotate`

CLI usage:

- `python setup_runtime.py --runtime-dir <...> --build --fetch-mc-annotate --json`
- `python setup_runtime.py --runtime-dir <...> --build --fetch-mc-annotate --fetch-mccore --json`

Important dependency note:

- MC-Annotate depends on `MCCORE` headers/libraries.
- Official MCCORE source is available at `https://github.com/major-lab/mccore`.
- `--fetch-mccore` lets setup clone/build/install MCCORE first, then build MC-Annotate against it.
- If MCCORE build fails on your platform/toolchain, setup reports the exact compiler/configure error.

This optional fetch/build step is setup-time only. Source-7 runtime behavior remains strict during runs.
