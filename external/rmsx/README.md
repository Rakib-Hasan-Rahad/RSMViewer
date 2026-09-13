# RNAMotifScanX Placeholder

Place the external RNAMotifScanX package here using this fixed layout:

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

Executables must have execute permission. RNAMotifScanX and its build/source tree are intentionally not distributed with RSMViewer.

## Preannotated data

Place the expanded preannotated dataset at:

```text
external/rmsx_preannotated/rmsx_work_default/<pdb_id>/
```

For example:

```text
external/rmsx_preannotated/rmsx_work_default/1s72/_prep_main/1S72.pdb
external/rmsx_preannotated/rmsx_work_default/1s72/0/1s72_0.rmsx.in
external/rmsx_preannotated/rmsx_work_default/1s72/0/k-turn_consensus.log
```

RSMViewer also accepts the packed fallback file `external/rmsx_preannotated/PDB_prebuild.tgz`. The expanded directory is preferred for cached target preparation.

Set behavior in `config/rmsx_config.json`:

```json
"data_mode": "preannotated"
```

Use `"run_from_scratch"` to execute the configured RMSX tools instead. Keep the research-paper cutoffs in `pvalue_thresholds`; P-values cannot be overridden from PyMOL commands.
