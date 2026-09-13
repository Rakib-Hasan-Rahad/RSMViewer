# RMSX Preannotated Data

RSMViewer loads RNAMotifScanX (Source 7) results from preannotated data placed
in this directory (`external/rmsx_preannotated/`).

## Recommended: paste the extracted folder (fast)

Extract `rmsx_preannotated_input_output.tar.gz` here so that the extracted
folder sits directly in this directory:

```
external/rmsx_preannotated/rmsx_work_default/
```

The extracted folder is named **`rmsx_work_default`** and contains one
subfolder per PDB (for example `rmsx_work_default/1s72/...`). When this folder
is present, RSMViewer reads from it directly, which is much faster than reading
the compressed archive.

To extract from this directory:

```bash
tar -xzf rmsx_preannotated_input_output.tar.gz
```

## Fallback: keep the compressed archive

If you forget to extract the folder, keep the compressed archive here with the
exact name:

```
external/rmsx_preannotated/rmsx_preannotated_input_output.tar.gz
```

RSMViewer automatically falls back to reading this archive when
`rmsx_work_default/` is missing. The first load extracts each PDB's results
once and caches them, so subsequent loads (and later PyMOL sessions) are fast.

## Notes

- RSMViewer always prefers the extracted `rmsx_work_default/` folder and only
  falls back to the `.tar.gz` archive when the folder is absent or has no data
  for the requested PDB.
- Leave this directory unchanged when preannotated data is not used.
- To regenerate results from scratch instead of using preannotated data, set
  `data_mode` to `run_from_scratch` in `config/rmsx_config.json`.
