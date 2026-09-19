# RNAMotifScanX (RMSX) Setup Guide

How to use RNAMotifScanX annotations in RSMViewer as a first-time user, on
**macOS**, **Windows** and **Linux**.

RSMViewer has two RMSX **data modes**, chosen with `data_mode` in
`config/rmsx_config.json`:

| `data_mode` | What it gives you | Extra software? | Internet? |
| --- | --- | --- | --- |
| `preannotated` (default) | Precomputed RMSX annotations for the PDB. | None | First use of each PDB |
| `run_from_scratch` | A fresh RNAMotifScanX run for the PDB. | Installed automatically on first use (see section 3) | First run only |

In both modes you load the annotations the same way:

```text
rmv_fetch 1S72
rmv_db RNAMotifScanX
```

Two housekeeping commands exist for RMSX:

| Command | Purpose |
| --- | --- |
| `rmv_setup RNAMotifScanX` | Optional. Prepares everything for `run_from_scratch` ahead of time (otherwise this happens automatically the first time you run). |
| `rmv_rmsx_doctor` | Health check: shows what is ready, what is missing, and how to fix it. |

---

## 1. First-time quick start

### Option A: precomputed annotations (nothing to install)

1. Keep `"data_mode": "preannotated"` in `config/rmsx_config.json` (the default).
2. In PyMOL:
   ```text
   rmv_fetch 1S72
   rmv_db RNAMotifScanX
   ```

The first time you request a PDB its annotations are downloaded automatically;
after that they load instantly and offline.

### Option B: run RNAMotifScanX yourself

1. Set `"data_mode": "run_from_scratch"` in `config/rmsx_config.json`.
2. Make sure your computer meets the requirements in section 3 (one-time, and
   only for macOS and Windows).
3. In PyMOL:
   ```text
   rmv_fetch 1S72
   rmv_db RNAMotifScanX
   ```

**The first time you run it, RSMViewer downloads all the required files and
third-party software by itself** from
[https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx.tar.gz](https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx.tar.gz)
into `external/rmsx/` (about 30 MB, checked against its `.sha256` checksum), then
prepares the scanner for your computer and runs it. This takes a few minutes and
prints its progress in the PyMOL window. Later runs start immediately. PyMOL
waits until the run finishes; 1S72 takes about a minute and a half.

If you would rather do the preparation first, run `rmv_setup RNAMotifScanX`, and
use `rmv_rmsx_doctor` at any time to see whether everything is ready.

---

## 2. Preannotated mode

```json
"data_mode": "preannotated"
```

To provide the most up-to-date RNAMotifScanX annotations, RSMViewer collects them
live from our server, **one PDB at a time**, the first time you request that PDB.

Where it looks, in order:

1. the local folder `external/rmsx_preannotated/rmsx_work_default/<pdb>/`;
2. the download link (the PDB ID is **lower case**):

```text
https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx_work_default/<pdb>.tar.gz

1S72  ->  https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx_work_default/1s72.tar.gz
```

A `1s72` sample is included, so `rmv_db RNAMotifScanX` works for 1S72 straight
away. You can also download an archive yourself (browser or `curl -O <link>`) and
extract it into `external/rmsx_preannotated/rmsx_work_default/`.

### Download criteria

A download is attempted only when all of these hold:

- `data_mode` is `preannotated`;
- the PDB ID is 4 characters long;
- there is no local `rmsx_work_default/<pdb>/` folder with results for it;
- the computer can reach the server (no proxy or firewall block).

The archive is accepted only if it extracts safely (no absolute or `..` paths, no
links or device files) and contains at least one `*_consensus.log` result file. It
is unpacked in a temporary folder and moved into place when complete, so an
interrupted download never leaves a broken folder. An existing folder is merged
into and nothing in it is deleted.

If the server has no archive for the PDB (HTTP 404) the dataset does not cover it
yet, and RSMViewer tells you so; use `run_from_scratch` for that structure. If the
download fails for another reason (offline, proxy, TLS), the message shows the
link so you can download it by hand.

### Getting newer data

A PDB that is already local is not downloaded again. To get the server's newest
version, delete its folder (`external/rmsx_preannotated/rmsx_work_default/<pdb>/`)
and run `rmv_db RNAMotifScanX` again. `rmv_reset cache` does not delete it.

---

## 3. Running RNAMotifScanX yourself (`run_from_scratch`)

```json
"data_mode": "run_from_scratch"
```

### 3.1 What you need first

| Your computer | One-time requirements |
| --- | --- |
| **Linux (x86-64)** | Nothing. |
| **Linux (ARM)** | `g++`, Boost and zlib: `sudo apt-get install g++ libboost-all-dev zlib1g-dev` (Debian/Ubuntu) or `sudo dnf install gcc-c++ boost-devel zlib-devel` (Fedora), **or** Docker. |
| **macOS (Apple Silicon or Intel)** | Xcode command line tools: `xcode-select --install`. [Homebrew](https://brew.sh) is used to install Boost automatically. Alternative with no compiler: `brew install docker colima`. |
| **Windows** | **WSL2** (recommended): in an administrator PowerShell run `wsl --install`, reboot, and finish the Linux user setup. Alternative: **Docker Desktop**. |
| **All** | An internet connection for the first run. |

macOS and Windows need a little more than Linux because the scanner is a Linux
program that is prepared for your computer on first use. RSMViewer picks the
right method for you (native build, WSL2 or Docker) and checks that the scanner
starts before using it.

An Intel build of PyMOL running under Rosetta on an Apple-silicon Mac is handled
automatically.

> Windows (WSL2 or Docker Desktop) has not been verified on a real Windows
> machine yet. If something fails there, please report it together with the
> output of `rmv_rmsx_doctor`.

### 3.2 Step by step (first-timer)

1. Install the requirements for your computer from the table above.
2. Open `config/rmsx_config.json` and set `"data_mode": "run_from_scratch"`.
3. Start PyMOL and run:
   ```text
   rmv_fetch 1S72
   rmv_db RNAMotifScanX
   ```
4. Wait. The first run downloads and prepares everything, then runs the
   analysis. Progress is printed as it goes.
5. Use the results as usual: `rmv_list`, `rmv_select`, `rmv_view`, and so on.

It also works in a combined load, e.g. `rmv_db RNA3DMotifAtlas,RNAMotifScanX`.

To scan only some chains, set `"scan_chains": ["9"]` in the config.

### 3.3 Check or prepare ahead of time

```text
rmv_rmsx_doctor          # what is ready / missing
rmv_setup RNAMotifScanX  # prepare everything now instead of at first run
```

`rmv_rmsx_doctor` ends with READY or NOT READY and, when something is missing,
the exact command or install step that fixes it.

---

## 4. Configuration reference (`config/rmsx_config.json`)

Relative paths are resolved against the folder containing the config file. The
file is re-read on every `rmv_db`, so edits take effect without restarting PyMOL.

| Key | Default | Purpose |
| --- | --- | --- |
| `data_mode` | `preannotated` | `preannotated` or `run_from_scratch`. |
| `scan_runtime` | `auto` | `auto` (native, then WSL2, then Docker), or force `native`, `wsl` or `docker`. |
| `rmsx_executable` | empty | Optional path to a `scan` program of your own; empty = the one RSMViewer prepares. |
| `scan_chains` | `[]` | Chains to run in `run_from_scratch`; empty = every chain of the PDB. |
| `pdb_prebuild_dir` | `../external/rmsx_preannotated/rmsx_work_default` | Local per-PDB data folder. |
| `preannotated_base_url` | `https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx_work_default` | Server folder with one `<pdb>.tar.gz` per PDB. |
| `rmsx_bundle_url` | `https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx.tar.gz` | Where the RNAMotifScanX files are downloaded from on first `run_from_scratch` use. |
| `query_motifs_dir` | empty | Directory of your own `<family>_consensus.struct` query files; empty = the standard set. |
| `output_dir` | `../output/rmsx_results` | Base folder for working results. |
| `motif_families` | all five | Families to load: `k-turn`, `c-loop`, `sarcin-ricin`, `reverse-kturn`, `e-loop`. |
| `num_threads` | `4` | Threads used by the scanner. |
| `scan_timeout_seconds` | `21600` | Time limit per scan (6 h). |
| `pvalue_thresholds` | see file | Per-family acceptance cutoff applied when annotations are loaded (both modes). |

Default P-value cutoffs: `KINK-TURN` 0.066, `C-LOOP` 0.044, `SARCIN-RICIN` 0.040,
`REVERSE-KINK-TURN` 0.018, `E-LOOP` 0.018 (0.05 for any other name).

Example: precomputed annotations for K-turns and sarcin-ricin only.

```json
{
  "data_mode": "preannotated",
  "motif_families": ["k-turn", "sarcin-ricin"]
}
```

Example: run it yourself with Docker, chain 9 only, 8 threads.

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

**First run cannot download the RNAMotifScanX files.** You may be offline or
behind a proxy. Download
[rmsx.tar.gz](https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx.tar.gz),
extract it, and copy its contents into `external/rmsx/` so that
`external/rmsx/RNAMotifScanX_src/` exists; then run again.

**"checksum mismatch".** The download was interrupted or altered. Run again; if it
repeats, download the archive manually as above.

**macOS: Boost or a compiler is missing.** Run `xcode-select --install` and
`brew install boost`, then `rmv_setup RNAMotifScanX` again, or use Docker:
`brew install docker colima`.

**Docker is installed but "daemon not reachable".** Start it (`colima start` on
macOS, Docker Desktop on Windows), or run `rmv_setup RNAMotifScanX`, which starts
Colima for you.

**Windows: "WSL is not installed".** In an administrator PowerShell run
`wsl --install`, reboot, then run `rmv_setup RNAMotifScanX`.

**"No preannotated RNAMotifScanX results are published for <PDB>" (HTTP 404).**
The server does not have that PDB yet. Use `run_from_scratch` for it.

**Preannotated download fails with a network or TLS error.** Download
`<preannotated_base_url>/<pdb>.tar.gz` (the link is printed in the message) with a
browser or `curl -O`, and extract it into
`external/rmsx_preannotated/rmsx_work_default/`.

**No accepted motifs.** Every P-value may exceed its cutoff (see
`pvalue_thresholds`). Zero hits is a valid result, not an error.

---

## 6. Reporting a result

Record: your operating system and CPU, `data_mode`, the date the RNAMotifScanX
files (or the PDB's archive) were downloaded, and the P-value cutoffs in
`config/rmsx_config.json`.
