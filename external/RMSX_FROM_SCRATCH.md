# RNAMotifScanX: Default and From-Scratch Workflows

This guide describes how RSMViewer uses RNAMotifScanX (RMSX), the role of `config/rmsx_config.json`, and how to run RMSX for the first time on a new computer.

## Default behavior: preannotated results

RSMViewer uses preannotated RMSX results by default:

```json
"data_mode": "preannotated"
```

With this setting, RSMViewer does not execute MC-Annotate, RNAVIEW, or the RMSX `scan` executable. It reads consensus logs from:

```text
external/rmsx_preannotated/rmsx_work_default/<pdb_id_lowercase>/
```

Example:

```text
external/rmsx_preannotated/rmsx_work_default/1s72/
external/rmsx_preannotated/rmsx_work_default/1kxk/
```

Each PDB directory may contain `_prep_main/` and one directory per RNA chain. Chain directories contain family logs such as:

```text
k-turn_consensus.log
c-loop_consensus.log
sarcin-ricin_consensus.log
reverse-kturn_consensus.log
e-loop_consensus.log
```

RSMViewer combines matching logs from all chains for each family and places them in the normal output cache. This is the recommended workflow on macOS and Windows.

```text
rmv_fetch 1S72
rmv_db RNAMotifScanX
```

To support another PDB, place its preannotated directory under `rmsx_work_default`, or configure the supplied preannotated archive.

## What from-scratch means

A complete from-scratch run has three stages:

```text
PDB structure
  -> MC-Annotate and RNAVIEW
  -> .rmsx.in and .rmsx.nch target files
  -> RNAMotifScanX scan for each motif family
  -> result_0_100_withbs.log files
```

The `_prep_main` files supplied in `rmsx_work_default` are already prepared inputs. Using them with `data_mode: "run_from_scratch"` skips MC-Annotate/RNAVIEW preparation and reruns only the core `scan` stage. This is useful for testing the scanner, but it is not a complete annotation-preparation run.

## Platform requirements

The original release executables are Linux binaries:

- RNAMotifScanX `scan`: Linux x86-64
- RNAVIEW: Linux x86-64
- MC-Annotate: Linux 32-bit Intel

The easiest first-time setup is a Linux x86-64 computer or VM. Windows users should use WSL2 or a Linux VM. Apple-Silicon macOS users need an x86-64 Linux VM/container with emulation; native arm64 macOS cannot execute these binaries directly.

The repository includes an optional scan wrapper:

```text
external/rmsx/bin/scan_docker_x86_64.sh
```

On Apple Silicon, install Docker, Colima, QEMU guest support, and start an x86-64 VM:

```bash
brew install colima docker lima-additional-guestagents
colima start --arch x86_64 --cpu 4 --memory 4 --disk 20 --vm-type qemu
chmod +x external/rmsx/bin/scan_docker_x86_64.sh
```

This wrapper runs the unchanged upstream `scan` executable. Emulation is substantially slower than native Linux, especially for large structures such as 1S72.

## Runtime layout

Obtain the official RNAMotifScanX release package. Do not modify the RMSX software. The runtime should provide:

```text
external/rmsx/
  RNAMotifScanX_src/
    scan
    Queries/*_consensus.struct
  bin/
    scan                 native Linux executable or scan wrapper
    MC-Annotate           required for complete preparation
    rnaview               required for complete preparation
  RNAVIEW/
    BASEPARS/
```

On Linux, make native executables executable:

```bash
chmod +x external/rmsx/bin/scan
chmod +x external/rmsx/bin/MC-Annotate
chmod +x external/rmsx/bin/rnaview
```

## Configuration file

`config/rmsx_config.json` controls the pipeline:

| Key | Purpose |
| --- | --- |
| `data_mode` | `preannotated` by default; `run_from_scratch` runs the scanner. |
| `rmsx_executable` | Path to `scan` or the Docker/QEMU wrapper. |
| `mc_annotate_executable` | Path to MC-Annotate for complete preparation. |
| `rnaview_executable` | Path to RNAVIEW for complete preparation. |
| `rnaview_dir` | RNAVIEW directory containing `BASEPARS`. |
| `query_motifs_dir` | Directory containing the five `*_consensus.struct` query files. |
| `pdb_prebuild_dir` | Optional extracted directory containing prepared `.rmsx.in`/`.nch` files. |
| `pdb_prebuild_archive` | Optional archive containing prepared targets. |
| `output_dir` | Destination for generated family logs. |
| `motif_families` | Families to scan: k-turn, c-loop, sarcin-ricin, reverse-kturn, e-loop. |
| `pvalue_thresholds` | Family-specific acceptance cutoffs used when annotations are loaded. |

Use absolute paths when launching the runner from a directory other than the repository root. Relative paths are resolved from the current working directory.

## True from-scratch configuration

For a complete preparation and scan run, use compatible executables and either provide a local PDB/CIF or enable downloads:

```json
{
  "data_mode": "run_from_scratch",
  "rmsx_executable": "/absolute/path/to/scan",
  "mc_annotate_executable": "/absolute/path/to/MC-Annotate",
  "rnaview_executable": "/absolute/path/to/rnaview",
  "rnaview_dir": "/absolute/path/to/RNAVIEW",
  "query_motifs_dir": "/absolute/path/to/Queries",
  "pdb_prebuild_dir": "",
  "pdb_prebuild_archive": "",
  "output_dir": "/absolute/path/to/output/rmsx_results",
  "auto_download_cif": true,
  "auto_download_pdb": true,
  "num_threads": 4,
  "max_strands": 3,
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

Leaving both prebuild settings empty forces the runner to perform MC-Annotate and RNAVIEW preparation before scanning. If prepared inputs are available and you want to reuse them, keep `pdb_prebuild_dir` configured; the preparation stage will be skipped.

## Run from the terminal

From the repository root:

```bash
cd "/absolute/path/to/RSMViewer-PyMOL PLUGIN"
python3 rsmviewer/tools/rmsx_runner.py \
  --config config/rmsx_config.json \
  --pdb 1KXK \
  --fresh
```

For the larger example:

```bash
python3 rsmviewer/tools/rmsx_runner.py \
  --config config/rmsx_config.json \
  --pdb 1S72 \
  --fresh
```

Check existing results without running:

```bash
python3 rsmviewer/tools/rmsx_runner.py \
  --config config/rmsx_config.json \
  --pdb 1S72 \
  --check
```

Expected outputs:

```text
output/rmsx_results/
  k-turn_consensus/result_0_100_withbs.log
  c-loop_consensus/result_0_100_withbs.log
  sarcin-ricin_consensus/result_0_100_withbs.log
  reverse-kturn_consensus/result_0_100_withbs.log
  e-loop_consensus/result_0_100_withbs.log
```

A valid family can have only a header and no accepted hits. That means the scan ran but no match passed the configured threshold; it is not automatically a runtime failure.

## Run through PyMOL

After choosing the RMSX source:

```text
rmv_fetch 1KXK
rmv_db RNAMotifScanX
```

The plugin invokes the configured runner and loads either preannotated or freshly generated family logs.

## Troubleshooting

### No runnable executable is configured

Check the configured path:

```bash
ls -l external/rmsx/bin/scan
```

On Apple Silicon, configure `scan_docker_x86_64.sh` and start Colima/QEMU first.

### Exec format error

A Linux x86-64 executable is being run directly on macOS or Windows. Use native Linux, WSL2, a Linux VM, or the Docker/QEMU wrapper.

### MC-Annotate or RNAVIEW fails

Check the PDB input, executable paths, and RNAVIEW `BASEPARS`. If valid `_prep_main` files exist, use them to test the core scanner separately.

### No accepted motifs

Inspect the generated `result_0_100_withbs.log` files. Header-only output means no motif passed the family P-value threshold.

## Reproducibility checklist

Record the operating system and architecture, RMSX release version, MC-Annotate/RNAVIEW versions, PDB and chain set, configuration file, whether the run used preannotated data, prepared inputs, or complete preparation, and the generated family logs.
