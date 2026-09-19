# RMSX Preannotated Data

RSMViewer loads RNAMotifScanX (Source 7) results from the preannotated data in
`external/rmsx_preannotated/rmsx_work_default/`, one subfolder per PDB (for
example `rmsx_work_default/1s72/...`).

## Automatic per-PDB download

You do not need to put anything here. To provide the most up-to-date
RNAMotifScanX annotations, RSMViewer collects the preannotated data live from
our server. The first time you request a PDB with `rmv_db RNAMotifScanX`,
RSMViewer downloads that PDB's results and extracts them into
`rmsx_work_default/<pdb>/`:

```text
<preannotated_base_url>/<pdb_lowercase>.tar.gz
e.g. https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx_work_default/1s72.tar.gz
```

The server address is `preannotated_base_url` in `config/rmsx_config.json`.
You can also download an archive yourself (browser or `curl -O <link>`) and
extract it into `rmsx_work_default/`.

### Download criteria

A download is attempted only when the PDB ID is 4 characters long, there are no
local results for it in `rmsx_work_default/<pdb>/`, and the server is reachable.
It is accepted only if the archive extracts safely (no absolute or `..` paths, no
links or device files) and contains at least one `*_consensus.log`. Extraction is
staged in a temporary folder and moved into place when complete, so an
interrupted download never leaves a broken folder. An existing folder is merged
into and nothing in it is deleted.

If the server has no archive for the PDB (HTTP 404), the dataset does not cover
it yet: RSMViewer says so, and you can scan the structure yourself with
`data_mode` set to `run_from_scratch` (see [../rmsx_setup.md](../rmsx_setup.md)).

## Lookup order

For each requested PDB, the first of these that has data is used:

1. `rmsx_work_default/<pdb>/` (already downloaded, or placed by you);
2. the download above.

A PDB that is already local is not re-downloaded. To get the newest results for
it, delete its folder and run `rmv_db RNAMotifScanX` again.

## Notes

- Downloaded folders are ignored by git, except the sample PDBs the repository
  tracks (`1s72`, `1ffk`; see `.gitignore`).
- Prepared `.rmsx.in`/`.rmsx.nch` inputs for `run_from_scratch` live in the same
  per-PDB folders; a download merges into an existing folder and never deletes
  your files.
- `rmv_reset cache` does not delete this folder: it is data, not cache.
- Full setup guide: [../rmsx_setup.md](../rmsx_setup.md).
