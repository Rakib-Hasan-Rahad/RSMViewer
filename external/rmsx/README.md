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
server (or read from the optional offline bundle; see
`external/rmsx_preannotated/README.md`).

## Run RMSX from scratch

Only use this directory when `data_mode` is set to `run_from_scratch` in
`config/rmsx_config.json`. Install the external RNAMotifScanX runtime in this
layout:

```text
external/rmsx/
  bin/
    scan
    MC-Annotate
    rnaview
  queries/
    *_consensus.struct
  RNAVIEW/
    BASEPARS/
```

The executables must have execute permission. RSMViewer does not version the
external runtime or its generated files.

## Configuration

```json
"data_mode": "preannotated"
```

With `preannotated`, RSMViewer prefers the extracted local folder, then downloads
the requested PDB's results from `preannotated_base_url`, and only then falls back
to `external/rmsx_preannotated/rmsx_preannotated_input_output.tar.gz`
(downloadable from Figshare:
[https://doi.org/10.6084/m9.figshare.33826795](https://doi.org/10.6084/m9.figshare.33826795)).
Keep the research-paper cutoffs in `pvalue_thresholds`; PyMOL commands do not
