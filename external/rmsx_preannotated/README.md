# RMSX Preannotated Data

RSMViewer loads RNAMotifScanX (Source 7) results from the preannotated data in
`external/rmsx_preannotated/rmsx_work_default/`, one subfolder per PDB (for
example `rmsx_work_default/1s72/...`).

## Automatic download (default)

You do not need to put anything here. To provide the most up-to-date
RNAMotifScanX annotations, RSMViewer collects the preannotated data live from
our server. The first time you request a PDB with
`rmv_db RNAMotifScanX`, RSMViewer downloads that PDB's results from the public
results server and extracts them into `rmsx_work_default/<pdb>/`:

```
<preannotated_base_url>/<pdb_lowercase>.tar.gz
e.g. https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx_work_default/1s72.tar.gz
```

Later loads of the same PDB read the local folder and do not download again.
This needs an internet connection for the first load of each PDB only. The
server address is `preannotated_base_url` in `config/rmsx_config.json`.

The preannotated dataset may not cover every PDB yet. If the server has no
results for a structure, RSMViewer tells you; you can then scan it yourself with
`data_mode` set to `run_from_scratch` (see `external/RMSX_FROM_SCRATCH.md`).

## Lookup order

For each requested PDB, the first of these that has data is used:

1. `rmsx_work_default/<pdb>/` (already downloaded, or placed by you);
2. download from the results server (above);
3. the optional compressed bundle (below), as an offline fallback.

## Optional: offline bundle

The full preannotated bundle `rmsx_preannotated_input_output.tar.gz` can also be
downloaded once from Figshare:
**[https://doi.org/10.6084/m9.figshare.33826795](https://doi.org/10.6084/m9.figshare.33826795)**.

Either extract it here so `rmsx_work_default/` holds every PDB, or keep the
archive at this exact path:

```
external/rmsx_preannotated/rmsx_preannotated_input_output.tar.gz
```

The archive is only read when a PDB is neither in `rmsx_work_default/` nor
available from the server. Reading it means decompressing the whole archive, so
it is slow; the extracted folder is much faster.

## Notes

- Downloaded folders are ignored by git, except the sample PDBs the repository
  tracks (`1s72`, `1ffk`, `4v88`; see `.gitignore`).
- Prepared `.rmsx.in`/`.rmsx.nch` inputs for `run_from_scratch` live in the same
  per-PDB folders; a download merges into an existing folder and never deletes
  your files.
- To regenerate results from scratch instead of using preannotated data, set
  `data_mode` to `run_from_scratch` in `config/rmsx_config.json`.
