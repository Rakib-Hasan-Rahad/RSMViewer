# RNAMotifScanX: Preannotated and Run-From-Scratch Workflows

This guide describes how RSMViewer uses RNAMotifScanX (RMSX), the role of
`config/rmsx_config.json`, and how a new user on **any platform** (Linux,
Apple-Silicon macOS, Intel macOS, or Windows) can run real RMSX scans.

RSMViewer supports **two** data modes, set by `data_mode` in
`config/rmsx_config.json`:

| `data_mode` | What it does | Executes RMSX? |
| --- | --- | --- |
| `preannotated` (default) | Reads bundled/downloaded consensus logs. Recommended for quick use. | No |
| `run_from_scratch` | Runs the real `scan` on the `.rmsx.in`/`.nch` inputs **you provide**. | Yes |

In `run_from_scratch` you generate the RMSX input files yourself — externally,
with RNAVIEW and MC-Annotate — and place them under `pdb_prebuild_dir`. RSMViewer
then runs only the RNAMotifScanX `scan` binary on those inputs and loads the
fresh output. **RSMViewer never runs MC-Annotate or RNAVIEW itself**, and it
never falls back to preannotated data for a from-scratch run.

If you only want annotations quickly, use `preannotated`.

---

## 1. Default behavior: preannotated results

```json
"data_mode": "preannotated"
```

Download the preannotated bundle (`rmsx_preannotated_input_output.tar.gz`) from
Figshare: **<https://doi.org/10.6084/m9.figshare.33826795>**, then extract it
inside `external/rmsx_preannotated/` so the `rmsx_work_default/` folder sits
there. RSMViewer does not execute any RMSX program in this mode. It reads
consensus logs from:

```text
external/rmsx_preannotated/rmsx_work_default/<pdb_id_lowercase>/
```

Each PDB directory contains one directory per RNA chain (plus `_prep_main/`),
holding family logs:

```text
k-turn_consensus.log
c-loop_consensus.log
sarcin-ricin_consensus.log
reverse-kturn_consensus.log
e-loop_consensus.log
```

In PyMOL:

```text
rmv_fetch 1S72
rmv_db RNAMotifScanX
```

This is the recommended workflow on macOS and Windows when you do not need a
live scan. To support another PDB, place its preannotated directory under
`rmsx_work_default`, or configure the supplied preannotated archive.

---

## 2. Run-from-scratch workflow: `run_from_scratch` (real scan on your inputs)

`run_from_scratch` runs the **unchanged upstream `scan` executable** directly on
the RMSX input files you provide and loads the freshly generated results. It:

- uses the `.rmsx.in` + matching `.rmsx.nch` files you place under
  `pdb_prebuild_dir` (see 2.1);
- **does not run MC-Annotate or RNAVIEW** — you generate the inputs externally;
- writes each run to its **own timestamped output directory**;
- captures the exact command, exit status, stdout, and stderr per scan;
- loads only that run's fresh output — it never substitutes, supplements, or
  falls back to preannotated data;
- reports a failed or cancelled run explicitly instead of showing false success.

The command it runs mirrors the header recorded in the preannotated logs:

```text
scan <family>_consensus.struct <pdb>_<chain>.rmsx.in \
     --map_pdb=<pdb>_<chain>.rmsx.nch \
     --pvalue 1.0 --num_threads <N> --write_alignment
```

### 2.1 Generate and place the input files

Generate the RMSX target inputs for your structure externally with RNAVIEW and
MC-Annotate (the RNAMotifScanX distribution ships the scripts that build
`.rmsx.in`/`.rmsx.nch` from a coordinate file). Then place one directory per
scanned chain under `pdb_prebuild_dir`:

```text
<pdb_prebuild_dir>/<pdb_id_lowercase>/<chain>/<pdb>_<chain>.rmsx.in
<pdb_prebuild_dir>/<pdb_id_lowercase>/<chain>/<pdb>_<chain>.rmsx.nch
```

For example, with the default `pdb_prebuild_dir`
(`external/rmsx_preannotated/rmsx_work_default`):

```text
external/rmsx_preannotated/rmsx_work_default/1kxk/A/1kxk_A.rmsx.in
external/rmsx_preannotated/rmsx_work_default/1kxk/A/1kxk_A.rmsx.nch
```

The `_prep_main/` staging folder that the RMSX tooling produces alongside the
per-chain folders is intermediate and is not scanned directly. If no matching
`.rmsx.in`/`.rmsx.nch` pair is found for a structure, the run stops with a clear
error — RSMViewer does **not** generate the inputs for you.

### 2.2 Enable it

```json
"data_mode": "run_from_scratch"
```

You do **not** need to set `rmsx_executable` to a working binary manually — the
runner resolves it automatically (see section 4).

### 2.3 Run from PyMOL

```text
rmv_fetch 1KXK
rmv_rmsx run 1KXK
```

Options:

```text
rmv_rmsx run <PDB> [CHAINS] [compare]
rmv_rmsx cancel
```

- `CHAINS`: optional, comma/space separated (for example `0,9`). Omit to scan
  every prepared chain.
- `compare`: optional development check that prints fresh-vs-preannotated hit
  counts. It never alters the loaded results.
- `cancel`: stops an in-progress run and **terminates the running
  scanner/container**, not just the remaining queue.

The scan runs off the PyMOL GUI thread, so PyMOL stays responsive. Progress is
reported per chain and family, and results load automatically when the run
finishes successfully.

Example status messages:

```text
Using your prepared RMSX inputs; MC-Annotate and RNAVIEW are not run by RSMViewer.
Running RNAMotifScanX for PDB 1KXK, chain A.
Loaded annotations from the newly generated RMSX output.
Results saved to: .../output/rmsx_results/run_from_scratch/1KXK_<timestamp>
```

### 2.4 Run from the terminal (standalone)

```bash
python rsmviewer/tools/rmsx_runner.py \
  --config config/rmsx_config.json \
  --pdb 1KXK --scan-prepared

# limit to specific prepared chains and choose an output directory
python rsmviewer/tools/rmsx_runner.py \
  --config config/rmsx_config.json \
  --pdb 1S72 --scan-prepared --chains 0,9 \
  --out output/rmsx_results/run_from_scratch/1S72
```

Use whichever Python launches PyMOL/RSMViewer on your machine (`python`,
`python3`, or an absolute interpreter path). The `--scan-prepared` flag is the
standalone entry point for the `run_from_scratch` scan.

### 2.5 Output layout for a run

```text
output/rmsx_results/run_from_scratch/<PDB>_<timestamp>/
  <family>_consensus/result_0_100_withbs.log   aggregated per family (loaded)
  _runs/chain_<c>/<family>_consensus/
    command.txt        exact command used
    scan.stdout.log    raw scan output
    scan.stderr.log    raw scan errors
```

A successful scan with **zero hits** is valid: `result_0_100_withbs.log` exists
with header-only content. A **failed** scan (non-zero exit) is reported as an
error and its family log is not written or loaded.

---

## 3. Platform setup for a real scan

The upstream `scan` executable is a **Linux x86-64** binary. How you run it
depends on your platform.

### 3.1 Native x86-64 Linux (simplest, fastest)

No Docker required. The runner auto-detects the bundled ELF at
`external/rmsx/RNAMotifScanX_src/scan`. Just ensure it is executable:

```bash
chmod +x external/rmsx/RNAMotifScanX_src/scan
python rsmviewer/tools/rmsx_runner.py --config config/rmsx_config.json --pdb 1KXK --scan-prepared
```

### 3.2 Apple-Silicon macOS / Intel macOS / Windows (Docker + Colima)

macOS and Windows cannot execute the Linux binary directly. RSMViewer ships a
wrapper that runs the unchanged binary inside an x86-64 Linux container:

```text
external/rmsx/bin/scan_docker_x86_64.sh
```

Install Docker + Colima and start an **x86_64** VM:

```bash
brew install colima docker
chmod +x external/rmsx/bin/scan_docker_x86_64.sh

# small structures (fast):
colima start --arch x86_64 --cpu 4 --memory 4

# large structures (e.g. 1S72 chain 0, ~2900 nt) need more memory:
colima start --arch x86_64 --cpu 6 --memory 10
```

Important notes:

- Ensure the Docker/Colima binaries are on your `PATH`. With Homebrew they live
  in `/opt/homebrew/bin` (Apple Silicon) or `/usr/local/bin` (Intel). If a
  terminal reports "no Docker", it is usually a `PATH` problem, not a missing
  install:

  ```bash
  export PATH="/opt/homebrew/bin:$PATH"
  ```

- Distinguish a missing binary from a stopped daemon:

  ```bash
  command -v docker colima        # installed?
  colima status                   # running?
  docker info                     # daemon reachable?
  ```

  If Colima is stopped, start it with `colima start`. You do not need to
  recreate the VM.

- The runner selects the wrapper automatically when a native binary cannot run
  and Docker is available; you do not have to edit `rmsx_executable`.

- Emulation is much slower than native Linux. Very large structures may exhaust
  memory or crash the scanner under emulation even with a 10 GB VM. If a large
  chain fails under emulation, run that structure on a native x86-64 Linux
  machine (section 3.1).

Windows: use WSL2 (native Linux, section 3.1) or Docker Desktop with a Linux
x86-64 container equivalent to the wrapper.

---

## 4. How the executable is resolved

For `run_from_scratch`, the runner tries, in order:

1. a **native** `scan` that can run on the current host:
   the configured `rmsx_executable`, then
   `external/rmsx/RNAMotifScanX_src/scan`, then `external/rmsx/bin/scan`;
2. the **Docker wrapper** `external/rmsx/bin/scan_docker_x86_64.sh` when Docker
   is installed.

If none is runnable, the run stops with an explicit reason (for example, a
Linux binary on macOS with no Docker) and does **not** fall back to
preannotated data.

---

## 5. Configuration file

`config/rmsx_config.json`:

| Key | Purpose |
| --- | --- |
| `data_mode` | `preannotated` or `run_from_scratch`. |
| `rmsx_executable` | Optional explicit `scan` path. Leave as bundled default; the runner also finds the ELF and wrapper automatically. |
| `query_motifs_dir` | Directory with the five `*_consensus.struct` query files. If empty, the runner finds `external/rmsx/RNAMotifScanX_src/Queries`. |
| `pdb_prebuild_dir` | Directory holding your prepared per-chain `.rmsx.in`/`.nch` inputs (used by `run_from_scratch`). |
| `pdb_prebuild_archive` | Optional archive of prepared inputs / preannotated logs. |
| `output_dir` | Base destination for generated logs. |
| `motif_families` | Families to scan: k-turn, c-loop, sarcin-ricin, reverse-kturn, e-loop. |
| `num_threads` | Threads passed to `scan` (default 4). |
| `pvalue_thresholds` | Family-specific acceptance cutoffs applied when annotations are loaded. |

Relative paths in the config are resolved relative to the config file's
directory.

Your prepared inputs live under `pdb_prebuild_dir`, one directory per chain:

```text
external/rmsx_preannotated/rmsx_work_default/1kxk/A/1kxk_A.rmsx.in
external/rmsx_preannotated/rmsx_work_default/1kxk/A/1kxk_A.rmsx.nch
external/rmsx_preannotated/rmsx_work_default/1s72/0/1s72_0.rmsx.in
external/rmsx_preannotated/rmsx_work_default/1s72/9/1s72_9.rmsx.nch
```

---

## 6. Generating the input files externally

`run_from_scratch` does **not** regenerate annotations inside RSMViewer. You
produce the `.rmsx.in`/`.rmsx.nch` target files yourself and drop them under
`pdb_prebuild_dir`; RSMViewer then runs only the `scan` step:

```text
PDB / mmCIF structure   (you)
  -> RNAVIEW + MC-Annotate        (you, external tools)
  -> .rmsx.in and .rmsx.nch       (you, placed under pdb_prebuild_dir)
  -> RNAMotifScanX scan           (RSMViewer, run_from_scratch)
  -> result_0_100_withbs.log      (RSMViewer output, loaded)
```

The RNAMotifScanX distribution ships the annotation/preparation scripts that
build `.rmsx.in`/`.rmsx.nch` from a coordinate file using RNAVIEW and
MC-Annotate. Run those on your machine (native Linux is simplest), then place
the resulting per-chain pairs as shown in section 2.1. RSMViewer intentionally
does not call MC-Annotate or RNAVIEW, and it does not download structures for
this mode — you control exactly which inputs are scanned.

Minimal `run_from_scratch` config:

```json
{
  "data_mode": "run_from_scratch",
  "rmsx_executable": "../external/rmsx/bin/scan",
  "pdb_prebuild_dir": "../external/rmsx_preannotated/rmsx_work_default",
  "output_dir": "../output/rmsx_results",
  "num_threads": 4,
  "motif_families": ["k-turn", "c-loop", "sarcin-ricin", "reverse-kturn", "e-loop"],
  "pvalue_thresholds": {
    "KINK-TURN": 0.066,
    "C-LOOP": 0.044,
    "SARCIN-RICIN": 0.040,
    "REVERSE-KINK-TURN": 0.018,
    "E-LOOP": 0.018
  }
}
```

Terminal:

```bash
python rsmviewer/tools/rmsx_runner.py --config config/rmsx_config.json --pdb 1KXK --scan-prepared
python rsmviewer/tools/rmsx_runner.py --config config/rmsx_config.json --pdb 1S72 --check
```

---

## 7. Development check: fresh vs preannotated

When validating a new setup, compare the freshly generated results with the
bundled preannotated dataset using equivalent templates and parameters:

```text
rmv_rmsx run 1KXK compare
```

Expect close but not necessarily identical counts. The bundled query templates
in `external/rmsx/RNAMotifScanX_src/Queries/*.struct` differ slightly from the
templates used to generate the preannotated dataset (and thread count may
differ), which can produce small per-family hit-count differences. Investigate
differences rather than assuming exact equality. The preannotated data is only a
reference; it never replaces the fresh results.

---

## 8. Troubleshooting

### "No Docker" but Docker is installed
`PATH` problem. Add the Homebrew bin directory and check the daemon:

```bash
export PATH="/opt/homebrew/bin:$PATH"
colima status || colima start
docker info
```

### Exec format error
A Linux x86-64 binary was run directly on macOS/Windows. Use native Linux,
WSL2, or the Docker/Colima wrapper (section 3.2).

### Scan fails on a very large structure under emulation
Increase VM memory (`colima start --arch x86_64 --cpu 6 --memory 10`). If it
still crashes, run that structure on a native x86-64 Linux machine.

### Missing or unreadable prepared inputs (run_from_scratch)
The run stops with "no prepared RMSX inputs found" when a structure has no
matching `.rmsx.in`/`.rmsx.nch` pair under `pdb_prebuild_dir`. Generate the
inputs externally (RNAVIEW + MC-Annotate) and place them as shown in section 2.1.
RSMViewer never generates them for you and never substitutes preannotated data.

### No accepted motifs
Header-only `result_0_100_withbs.log` means the scan ran but no motif passed the
family P-value threshold. This is a valid result, not a failure.

---

## 9. Reproducibility checklist

Record: operating system and architecture; whether the run used native Linux or
the Docker/Colima wrapper; RMSX release version (and MC-Annotate/RNAVIEW for
`run_from_scratch`); PDB and chain set; the `data_mode`; the configuration file;
the exact command from `command.txt`; and the generated family logs.
