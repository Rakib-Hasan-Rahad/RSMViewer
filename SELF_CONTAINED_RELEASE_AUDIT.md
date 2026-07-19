# Self-Contained Release Audit

Release definition used in this audit:

A user can clone/download the GitHub repository, install the plugin into PyMOL, and use all local functionality without installing, cloning, linking, configuring, or manually locating any external software.

This audit does **not** treat intentionally online features as packaging blockers.

---

## Executive Summary

Under this release definition, the repository is **not yet fully self-contained**, but the blocker set is much smaller than a full offline-everything interpretation.

### Must Fix Before Release

1. The published repository does not yet include the vendored FR3D runtime source tree required by Source 5.
2. Only `macos-arm64` RMSX runtime binaries are currently bundled.
3. Release packaging still needs to publish the missing local runtime assets instead of keeping them only in the local workspace.

### Platform Support Gap

1. Missing `macos-x86_64` RMSX binaries.
2. Missing `linux-x86_64` RMSX binaries.
3. Missing `windows-x86_64` RMSX binaries.

### Online Feature

1. `rmv_fetch <PDB_ID>` downloads a structure from the internet.
2. Source 3 uses the BGSU API.
3. Source 4 uses the Rfam API.

These are not packaging problems.

### Licensing Review

1. FR3D redistribution should be documented cleanly in release notices.
2. RNAMotifScanX redistribution status should be verified explicitly.
3. MC-Annotate redistribution/license terms should be verified explicitly.
4. MCCORE redistribution/license compliance should be verified explicitly.

---

## 1. Must Fix Before Release

## 1.1 Vendored FR3D runtime source is complete locally but not yet published

Current state:

- Source 5 local runtime expects an FR3D source tree under:
  - `rsmviewer/tools/fr3d_runtime/python/fr3d-python/`
- that tree exists in the local workspace and was verified as a full upstream checkout
- local verification found required files present for runtime use
- that tree is not part of the currently published GitHub repository (`git ls-files` count is zero under this path)

Why this is a release blocker:

- Source 5 cannot be zero-setup if the user has to clone FR3D separately

What was implemented now:

- `rsmviewer/tools/fr3d_runtime/setup_runtime.py` was changed so release setup no longer clones FR3D from the network
- Source 5 now expects the bundled runtime to exist in the package

What still remains:

1. publish the cleaned FR3D source tree into the repository
2. exclude nested `.git`, `__pycache__`, `venv`, and packaging leftovers

## 1.2 FR3D runtime clone-on-first-run behavior was a blocker (fixed)

Previous state:

- setup script could clone FR3D using `git clone`

Current state after implementation:

- fixed in code
- setup script now refuses external source acquisition and expects the bundled runtime tree

Files changed:

- `rsmviewer/tools/fr3d_runtime/setup_runtime.py`
- `rsmviewer/gui.py`

## 1.3 Missing published local runtime assets

Current state:

- important local runtime assets still exist only in the local workspace and are not published yet:
  - `rsmviewer/tools/fr3d_runtime/python/fr3d-python/`
  - `rsmviewer/tools/rmsx_runtime/src/`

Why this matters:

- for a self-contained local-feature release, everything needed by the packaged local runtime must live inside the repository

Important note:

- `rsmviewer/tools/rmsx_runtime/src/` is not necessarily required for end-user runtime if all finished binaries are bundled
- it is mainly relevant if you want the repository to contain source provenance or future rebuild capability

---

## 2. Platform Support Gap

## 2.1 RMSX currently bundled only for macOS Apple Silicon

Bundled now:

- `rsmviewer/tools/rmsx_runtime/bin/macos-arm64/scan`
- `rsmviewer/tools/rmsx_runtime/bin/macos-arm64/MC-Annotate`
- `rsmviewer/tools/rmsx_runtime/bin/macos-arm64/libmccore.2.0.0.dylib`

Missing:

- `rsmviewer/tools/rmsx_runtime/bin/macos-x86_64/`
- `rsmviewer/tools/rmsx_runtime/bin/linux-x86_64/`
- `rsmviewer/tools/rmsx_runtime/bin/windows-x86_64/`

Why this is not a generic architecture bug but a packaging gap:

- runtime detection already exists in `rsmviewer/gui.py`
- the folder layout for automatic platform selection already exists
- the missing part is the platform-specific bundled executable set

What can be implemented only when assets exist:

1. publish the matching binaries into those folders
2. include dependent shared libraries where required
3. validate each platform with `rmv_rmsx_doctor`, `rmv_db 7`, and `rmv_load_motif`

---

## 3. Online Feature

These are intentionally online features and are **not** packaging blockers.

## 3.1 `rmv_fetch <PDB_ID>`

- downloads structures from the internet
- this is expected runtime behavior
- it does not mean the plugin is not self-contained from a packaging standpoint

## 3.2 Source 3

- uses the BGSU online source
- intentionally network dependent

## 3.3 Source 4

- uses the Rfam online source/API
- intentionally network dependent

Release implication:

- these should be documented as online features
- they should not be used to reject the local-runtime packaging goal

---

## 4. Licensing Review

These items are important for release quality and redistribution safety, but they are separate from the pure technical packaging problem.

## 4.1 FR3D

Evidence found:

- `rsmviewer/tools/fr3d_runtime/python/fr3d-python/pyproject.toml` declares `Apache-2.0`

Assessment:

- likely bundleable
- release should still include proper third-party notice/license documentation

## 4.2 RNAMotifScanX

Evidence found:

- source README exists
- inspected root metadata did not show a clear standalone license file

Assessment:

- redistribution permission should be verified explicitly before final release packaging

## 4.3 MC-Annotate

Evidence found:

- source/build metadata references GPL-style packaging
- source tree contains copyright notices

Assessment:

- redistribution terms should be explicitly verified and documented before broad binary bundling

## 4.4 MCCORE

Evidence found:

- source headers contain LGPL language
- some packaging metadata references GPLv3+

Assessment:

- license/compliance posture should be reviewed before final all-platform binary redistribution

---

## 5. What Was Implemented in This Re-Audit Pass

Implemented now:

1. Reclassified blockers using the requested release definition.
2. Removed FR3D external clone/bootstrap behavior from the setup path.
3. Updated FR3D user/runtime messages in `gui.py` to describe the bundled-runtime model.
4. Added `requirements-plugin.txt` to document standard FR3D Python dependencies.

Not yet implementable from the current published repo alone:

1. missing RMSX binaries for `macos-x86_64`, `linux-x86_64`, and `windows-x86_64`
2. publishing the cleaned vendored FR3D runtime source tree
3. publishing any source provenance/runtime source bundles you decide to keep for RMSX toolchain

---

## 6. Final Verdict

Using the refined release definition, the repository is close to the target model but is **not release-complete yet**.

The real remaining technical blockers are:

1. publish the bundled FR3D source runtime into the repository
2. publish all supported RMSX platform binaries into the runtime layout

The online features are not blockers.

The Python dependencies are acceptable as normal plugin dependencies if they are documented and installed through the standard plugin installation process.
