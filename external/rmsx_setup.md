# RNAMotifScanX (RMSX) Setup Guide

This guide explains how RSMViewer gets RNAMotifScanX annotations, what
`config/rmsx_config.json` controls, and how to set everything up on **macOS**,
**Windows** and **Linux**.

RSMViewer has two RMSX **data modes**. You choose one with `data_mode` in
`config/rmsx_config.json`:

| `data_mode` | What happens | Needs the scanner? | Needs internet? |
| --- | --- | --- | --- |
| `preannotated` (default) | Loads precomputed RMSX results for the PDB, downloaded per PDB from the project's server. | No | First use of each PDB |
| `run_from_scratch` | Runs the real `scan` program on the PDB's prepared inputs and loads the fresh output. | Yes | Only to fetch missing inputs |

Whichever mode you pick, the PyMOL commands are the same:

```text
rmv_fetch 1S72
rmv_db RNAMotifScanX
```

Only two commands exist for RMSX housekeeping:

| Command | Purpose |
| --- | --- |
| `rmv_setup RNAMotifScanX` | One-time setup of the scanner for `run_from_scratch` (not needed for `preannotated`). |
| `rmv_rmsx_doctor` | Read-only health check: which runtime would be used, what is missing, whether the server is reachable. |

Everything else (choosing the mode, threads, chains, P-value cutoffs) is done in
the config file (section 4).

> RSMViewer never runs MC-Annotate or RNAVIEW itself, and a failed
> `run_from_scratch` run is never silently replaced by preannotated data.

---

## 1. Quick start

**I just want to see RMSX annotations** (most users):

1. Leave `"data_mode": "preannotated"` in `config/rmsx_config.json`.
2. In PyMOL:
   ```text
   rmv_fetch 1S72
   rmv_db RNAMotifScanX
   ```
3. Nothing else to install. The first time you request a PDB, its results are
   downloaded automatically.

**I want to run RNAMotifScanX myself:**

1. Set `"data_mode": "run_from_scratch"` in `config/rmsx_config.json`.
2. In PyMOL: `rmv_setup RNAMotifScanX` (once), then `rmv_rmsx_doctor` to confirm
   it says READY.
3. `rmv_fetch 1S72` then `rmv_db RNAMotifScanX`.

Platform prerequisites for step 2 are in section 3.

---

## 2. Preannotated mode

```json
"data_mode": "preannotated"
```

No RMSX program runs. To provide the most up-to-date RNAMotifScanX annotations,
RSMViewer collects the precomputed results **live from our server, one PDB at a
time**, instead of shipping a large offline dataset.

### 2.1 Where the data comes from

For a requested PDB, RSMViewer looks in this order and uses the first source that
has data:

1. **Local folder**: `external/rmsx_preannotated/rmsx_work_default/<pdb>/`
2. **Download**: `<preannotated_base_url>/<pdb>.tar.gz`, extracted into the folder
   above and then read from there.

The download link has this form (the PDB ID is **lower case**):

```text
https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx_work_default/<pdb>.tar.gz

1S72  ->  https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx_work_default/1s72.tar.gz
```

You can also download an archive yourself (browser or `curl -O <link>`) and
extract it into `external/rmsx_preannotated/rmsx_work_default/`; it is then
treated as already downloaded.

### 2.2 Download criteria

A download is attempted only when **all** of these hold:

- `data_mode` is `preannotated` (or, in `run_from_scratch`, prepared inputs for
  the PDB are missing locally; see section 3.6);
- the PDB ID is a 4-character ID;
- there is no local `rmsx_work_default/<pdb>/` folder with results for it;
- the machine can reach the server (no proxy or firewall block).

A downloaded archive is accepted only if it extracts safely (no absolute or
`..` paths, no links or device files) and contains at least one
`*_consensus.log` result file. It is unpacked in a temporary folder and moved
into place only when complete, so an interrupted download never leaves a broken
half-folder behind. If a folder for that PDB already exists, the new files are
merged into it and nothing of yours is deleted.

If the server has **no archive for the PDB** (HTTP 404) the dataset does not
cover it yet: RSMViewer tells you so, and you can scan the structure yourself
with `run_from_scratch`. If the download fails for another reason (offline,
proxy, TLS), the error is printed with the link so you can download it by hand.

### 2.3 Keeping data current

A PDB that is already in the local folder is **not** re-downloaded. To get the
server's newest results for it, delete its folder
(`external/rmsx_preannotated/rmsx_work_default/<pdb>/`) and run
`rmv_db RNAMotifScanX` again. `rmv_reset cache` does not delete this folder.

### 2.4 What a PDB folder contains

One folder per RNA chain (plus `_prep_main/`) holding the family logs
(`k-turn_consensus.log`, `c-loop_consensus.log`, `sarcin-ricin_consensus.log`,
`reverse-kturn_consensus.log`, `e-loop_consensus.log`) and the prepared scanner
inputs (`<pdb>_<chain>.rmsx.in` and `.rmsx.nch`). Preannotated mode reads the
logs; `run_from_scratch` reads the inputs.

---

## 3. Run-from-scratch mode

```json
"data_mode": "run_from_scratch"
```

RSMViewer runs the real `scan` program for each chain and family of the PDB, then
loads the accepted hits exactly as it would preannotated ones.

### 3.1 Why setup is needed

RNAMotifScanX is a C++ program. The binary shipped in `external/rmsx/` is a
**Linux x86-64** executable, so it cannot run directly on macOS or Windows.
`rmv_setup RNAMotifScanX` picks the first of these that works on your machine
and verifies it by actually starting the scanner:

| Runtime | Platform | What setup does |
| --- | --- | --- |
| **native** | Linux x86-64 | Uses the bundled binary as is. |
| **native** | macOS, other Linux | **Builds `scan` from the bundled source** into `external/rmsx/bin/<platform>/` (C++ compiler + Boost). On macOS it runs `brew install boost` if Boost is missing. No virtualisation. |
| **WSL2** | Windows | Copies the Linux binary into your WSL distribution and runs it there. |
| **Docker** | any | Runs the Linux binary in an `ubuntu:22.04` `linux/amd64` container. On macOS, if Colima is installed but stopped, setup starts it. |

`scan_runtime: "auto"` (the default) tries them in the order native, WSL2,
Docker. You can force one (section 4).

### 3.2 macOS (Apple Silicon or Intel)

Recommended route: build natively (fastest, nothing runs in a VM).

```bash
xcode-select --install          # C++ compiler (skip if already installed)
brew install boost              # setup runs this for you if Homebrew exists
```

Then in PyMOL:

```text
rmv_setup RNAMotifScanX
rmv_rmsx_doctor
```

Alternative (no compiler): Docker through Colima.

```bash
brew install docker colima
colima start
```

then `rmv_setup RNAMotifScanX` (it starts Colima for you if it is stopped).

### 3.3 Linux

- **x86-64 (Intel/AMD):** nothing to install. The bundled binary is used as is
  (it needs only glibc). Run `rmv_setup RNAMotifScanX` once to verify it.
- **ARM Linux, or a binary that will not start:** install a compiler, Boost and
  zlib and let setup build the scanner:
  ```bash
  sudo apt-get install g++ libboost-all-dev zlib1g-dev     # Debian/Ubuntu
  sudo dnf install gcc-c++ boost-devel zlib-devel          # Fedora/RHEL
  ```
  or use Docker (`sudo apt-get install docker.io`, and add your user to the
  `docker` group).

### 3.4 Windows

The scanner cannot be built natively on Windows. Use **one** of:

- **WSL2** (recommended): in an *administrator* PowerShell run
  `wsl --install`, reboot, and finish the Linux user creation. Then in PyMOL run
  `rmv_setup RNAMotifScanX`; setup copies the scanner into your WSL distribution.
- **Docker Desktop**: install it, start it, and make sure it is running Linux
  containers. Then `rmv_setup RNAMotifScanX`.

> The Windows routes are implemented and unit-tested but have not been verified on
> a real Windows machine. If something fails, please report it together with the
> output of `rmv_rmsx_doctor`.

### 3.5 Verify and run

```text
rmv_rmsx_doctor
```

The doctor reports the platform, which runtime would be used, the build
toolchain, the scoring matrices and query models, prepared inputs for the loaded
PDB, and whether the results server is reachable. It ends with READY or NOT
READY and, when not ready, the exact command that fixes each problem.

Then:

```text
rmv_fetch 1S72
rmv_db RNAMotifScanX
```

It also works in a combined load, e.g. `rmv_db RNA3DMotifAtlas,RNAMotifScanX`.

**PyMOL pauses until the scan finishes** (as it does for FR3D), printing progress
per chain and family. On an Apple M4, all five families on both chains of 1S72
take about 80 seconds; larger structures take longer. A repeated `rmv_db` in the
same session reuses the scan; `rmv_refresh` or `rmv_reset session` forces a new
one.

### 3.6 Where the inputs come from

The scanner reads prepared inputs from `pdb_prebuild_dir`
(`external/rmsx_preannotated/rmsx_work_default`):

```text
<pdb_prebuild_dir>/<pdb>/<chain>/<pdb>_<chain>.rmsx.in
<pdb_prebuild_dir>/<pdb>/<chain>/<pdb>_<chain>.rmsx.nch
```

If the PDB has no prepared inputs locally, RSMViewer downloads
`<pdb>.tar.gz` (section 2.1, same link and criteria) **only to obtain these
inputs**; the precomputed logs inside are never used as scan results. For a
structure the server does not have, generate its `.rmsx.in`/`.rmsx.nch` files
yourself with RNAVIEW and MC-Annotate (the RNAMotifScanX distribution ships the
scripts) and place them at the paths above.

To scan only some chains set `"scan_chains": ["9"]` in the config.

### 3.7 What is run

For each chain and family, on the chosen runtime:

```text
scan <family>_consensus.struct <pdb>_<chain>.rmsx.in \
     --map_pdb=<pdb>_<chain>.rmsx.nch \
     --pvalue 1.0 --num_threads <N> --write_alignment
```

with `RNAMOTIFSCANX_PATH` pointing at `external/rmsx/RNAMotifScanX_src` (scoring
matrices). File paths are translated automatically for WSL2 and Docker. Query
models come from `Queries/reduced/` first: those reproduce the published
preannotated results (checked against the 1S72 logs), while the full `Queries/`
set also reports additional marginal hits.

### 3.8 Output

Each run writes to its own timestamped folder under `output_dir`:

```text
output/rmsx_results/run_from_scratch/<PDB>_<YYYYMMDD_HHMMSS>/
├── <family>_consensus/result_0_100_withbs.log   aggregated over chains (loaded)
└── _runs/chain_<c>/<family>_consensus/
    ├── command.txt        the exact command
    ├── scan.stdout.log    raw scanner output
    └── scan.stderr.log
```

### 3.9 Things to know about RNAMotifScanX itself

- **P-values are random estimates.** `scan` estimates each P-value by a random
  simulation seeded with the clock. Scores and alignments repeat exactly, but
  P-values (and so which borderline hits pass a cutoff such as 0.018) differ
  slightly from run to run and from the published preannotated data.
- **The scanner can crash on large RNAs.** The E-loop scan of the 2,900-nt 23S
  rRNA chain of 1S72 segfaults; the published preannotated log ends with
  `# scan failed rc=-11` too. RSMViewer keeps every complete alignment found
  before the crash, drops the truncated last one, ends the log with the same
  `# scan failed` marker, and warns you that later hits are missing. A scan that
  crashes with no complete alignment is reported as failed and nothing is loaded.

---

## 4. Configuration reference (`config/rmsx_config.json`)

Relative paths are resolved against the folder containing the config file. The
file is re-read on every `rmv_db`, so edits take effect without restarting PyMOL.

| Key | Default | Purpose |
| --- | --- | --- |
| `data_mode` | `preannotated` | `preannotated` or `run_from_scratch` (sections 2 and 3). |
| `scan_runtime` | `auto` | `auto` (native, then WSL2, then Docker), or force `native`, `wsl` or `docker`. |
| `rmsx_executable` | empty | Optional path to a `scan` binary of your own; empty = bundled or built one. |
| `scan_chains` | `[]` | Chains to scan in `run_from_scratch`; empty = every prepared chain. |
| `pdb_prebuild_dir` | `../external/rmsx_preannotated/rmsx_work_default` | Where downloaded results and prepared inputs live. |
| `preannotated_base_url` | `https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx_work_default` | Server folder holding `<pdb>.tar.gz` per PDB. |
| `query_motifs_dir` | empty | Directory of `<family>_consensus.struct` query files; empty = `external/rmsx/RNAMotifScanX_src/Queries/reduced`. |
| `output_dir` | `../output/rmsx_results` | Base folder for working results. |
| `motif_families` | all five | Families to load/scan: `k-turn`, `c-loop`, `sarcin-ricin`, `reverse-kturn`, `e-loop`. |
| `num_threads` | `4` | Threads passed to `scan`. |
| `scan_timeout_seconds` | `21600` | Per-scan time limit (6 h). |
| `pvalue_thresholds` | see file | Per-family acceptance cutoff applied when annotations are loaded (both modes). |

Default P-value cutoffs: `KINK-TURN` 0.066, `C-LOOP` 0.044, `SARCIN-RICIN` 0.040,
`REVERSE-KINK-TURN` 0.018, `E-LOOP` 0.018 (fallback 0.05 for other names).

Example: preannotated, K-turns and sarcin-ricin only.

```json
{
  "data_mode": "preannotated",
  "motif_families": ["k-turn", "sarcin-ricin"]
}
```

Example: run from scratch on Docker, only chain 9, 8 threads.

```json
{
  "data_mode": "run_from_scratch",
  "scan_runtime": "docker",
  "scan_chains": ["9"],
  "num_threads": 8
}
```

---

## 5. Troubleshooting

Start with `rmv_rmsx_doctor`; it names the exact problem.

**"no RNAMotifScanX scanner is available on this machine".** Run
`rmv_setup RNAMotifScanX`. The message lists why each runtime was rejected.

**macOS: setup says Boost or a compiler is missing.** Run
`xcode-select --install` and `brew install boost`, then `rmv_setup RNAMotifScanX`
again, or use Docker: `brew install docker colima`.

**macOS: "Platform: Darwin x86_64" on an Apple-silicon Mac, or "linking failed".**
Your PyMOL is an Intel build running under Rosetta. This is handled: RSMViewer
detects the real hardware (`rmv_rmsx_doctor` shows `arm64` and a Rosetta note),
builds the scanner natively for arm64 into `external/rmsx/bin/macos-arm64/`
(matching Homebrew's arm64 Boost), and runs it from PyMOL. If you still see the
error, update RSMViewer and rerun `rmv_setup RNAMotifScanX`; the compile and link
log is in `external/rmsx/bin/macos-arm64/build.log`.

**Docker is installed but "daemon not reachable".** Start it (`colima start` on
macOS, Docker Desktop on Windows) or run `rmv_setup RNAMotifScanX`, which starts
Colima for you.

**Windows: "WSL is not installed".** In an administrator PowerShell run
`wsl --install`, reboot, then `rmv_setup RNAMotifScanX`.

**"No preannotated RNAMotifScanX results are published for <PDB>" (HTTP 404).**
The server does not have that PDB yet. Use `run_from_scratch` for it.

**Download fails with a network or TLS error.** Download
`<preannotated_base_url>/<pdb>.tar.gz` (link printed in the message) with a
browser or `curl -O`, and extract it into
`external/rmsx_preannotated/rmsx_work_default/`.

**"could not download prepared inputs" in `run_from_scratch`.** The server has no
inputs for the PDB, or you are offline. Place the `.rmsx.in`/`.rmsx.nch` files
under `pdb_prebuild_dir` yourself (section 3.6).

**No accepted motifs after a successful scan.** Every P-value may exceed its
cutoff (see `pvalue_thresholds`). Zero hits is a valid result, not a failure.

---

## 6. Reproducibility checklist

When reporting an RMSX result, record: operating system and CPU; `data_mode`; for
`run_from_scratch` the runtime used (shown by `rmv_rmsx_doctor` and at the start
of each run), the RMSX release and the query directory (`reduced` or full); for
`preannotated` the date the PDB's archive was downloaded; and the P-value cutoffs
in `config/rmsx_config.json`.
