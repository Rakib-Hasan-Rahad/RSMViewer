# RNAMotifScanX Runtime

RSMViewer uses preannotated RNAMotifScanX results by default. In that mode it
does not execute the programs in this directory. The preannotated data is read
from `external/rmsx_preannotated/`. To provide the most up-to-date annotations it
is collected live from our server, per PDB, the first time it is requested; see
its README for details.

## Included sample

The repository includes a small extracted `1s72` sample at:

```text
external/rmsx_preannotated/rmsx_work_default/1s72/
```

It contains the consensus result logs for chains `0` and `9` and lets the
default `rmv_db RNA3DMotifAtlas,RNAMotifScanX` workflow run immediately for
`1S72`. Other structures are downloaded on first use from the public results
server, one PDB at a time (see `external/rmsx_preannotated/README.md`).

## Run RMSX from scratch

Only needed when `data_mode` is `"run_from_scratch"` in `config/rmsx_config.json`.
This directory holds everything required:

```text
external/rmsx/
  RNAMotifScanX_src/     C++ source, scoring matrices (mat/), query models
                         (Queries/reduced/, Queries/), and a Linux x86-64 build
                         of `scan`
  bin/<platform>/scan    a `scan` built from that source on your machine
                         (created by `rmv_setup RNAMotifScanX`; git-ignored)
```

The bundled `scan` is a Linux x86-64 executable. Run `rmv_setup RNAMotifScanX`
once: it uses that binary on Linux x86-64, builds one from source on macOS and
other Linux (needs a C++ compiler and Boost), or falls back to WSL2 (Windows) or
Docker. `rmv_rmsx_doctor` shows what is available. See
[../rmsx_setup.md](../rmsx_setup.md).

## Configuration

```json
"data_mode": "preannotated"
```

With `preannotated`, RSMViewer uses the local folder
`external/rmsx_preannotated/rmsx_work_default/<pdb>/` if it has results, and
otherwise downloads that PDB's `<pdb>.tar.gz` from `preannotated_base_url`.
Keep the research-paper cutoffs in `pvalue_thresholds`; PyMOL commands do not
