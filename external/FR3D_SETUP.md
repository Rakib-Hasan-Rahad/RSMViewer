# Running FR3D through RSMViewer

A short guide to set up and run the official BGSU **fr3d-python** from the
RSMViewer PyMOL plugin on any PC (Windows / macOS / Linux).

---

## 1. One-time setup

RSMViewer supports two FR3D modes. The repository default is `cache`, which
loads prepared local data from `external/fr3d/fr3d_cache/` and does not run
FR3D or make a network call. To execute the official FR3D Python pipeline,
set this field in `config/fr3d_config.json`:

```json
"data_mode": "run_from_scratch"
```

After changing the mode, run `rmv_db FR3D` again.

1. **Paste the FR3D code.** Put the official `fr3d-python` checkout here:

   ```
   external/fr3d/fr3d-python-latest/
   ```

   It must contain `fr3d/__init__.py` and `fr3d/search/FR3D.py`.

2. **Set up and check FR3D** from inside PyMOL when using
  `run_from_scratch` (uses PyMOL's own Python — no
  separate Python needed):

   ```
   rmv_setup FR3D
   ```

  This installs `numpy`, `scipy`, and `mmcif-pdbx` when needed, registers
  FR3D, then prints one concise status report. It says **ready** when the
  checkout, Python, dependencies, queries, and output folder are usable; if
  not, it names the missing requirement and the next command to run.

---

## 2. Run FR3D

```
rmv_fetch 1S72
rmv_db FR3D
```

With `data_mode: "cache"`, `rmv_db FR3D` loads the matching local cache file.
With `data_mode: "run_from_scratch"`, it annotates the loaded structure by
running the motif queries listed in the config. Results load straight into
RSMViewer in both modes.

> **First run is slow, later runs are fast.** FR3D's `geometric_*` queries
> embed large reference structures (e.g. `4V9F`, `7K00`, `8GLP`) in their
> definitions. FR3D must download and fully annotate those the first time,
> which can take several minutes and makes PyMOL appear frozen. The result is
> saved to a **persistent cache** that is reused by every later run, so
> re-running FR3D — on the same or a different structure — is fast.
>
> The cache lives at `output/fr3d_runs/_fr3d_cache/`. Deleting it only forces
> the one-time rebuild again; it is safe to keep.

Good structures to try: `1S72`, `4V9F`, `4V88`, `1HR2`, `1KXK`, `3CC2`,
`1FFK`, `1NBS`, `1JJ2`, `1Y0Q`, `2GIS`.

---

## 3. The config file — `config/fr3d_config.json`

```json
{
  "data_mode": "cache",
  "fr3d_python_path": "../external/fr3d/fr3d-python-latest",
  "cache_path": "../external/fr3d/fr3d_cache",
  "query_path": "../external/fr3d/fr3d-python-latest/fr3d/search/queries",
  "query_selection": "selected",
  "query_families": [
    "geometric_5_sarcin_ricin",
    "geometric_5_kink_turn_65553",
    "symbolic_6_GNRA_hairpin_symbolic",
    "symbolic_4_internal_loop",
    "symbolic_8_GU_tandem_at_end_of_helix"
  ],
  "run_output_path": "../output/fr3d_runs",
  "allow_network": true
}
```

| Key | What it does |
| --- | --- |
| `fr3d_python_path` | Where you pasted the FR3D checkout. |
| `query_path` | FR3D's own query folder (the motif definitions). |
| `query_selection` | Which queries to run (see table below). |
| `query_families` | The exact query file names to run when `query_selection` is `selected`. |
| `allow_network` | `true` lets FR3D fetch a structure / a query's reference template when needed. |
| `run_output_path` | Where results and provenance are written. |

### `query_selection` options

| Value | Runs |
| --- | --- |
| `selected` | **Default.** Only the queries listed in `query_families`. |
| `all` | Every query in `query_path` — the **full pipeline** (slow; many queries). |

---

## 4. Choose your own motifs

The motif queries live in:

```
external/fr3d/fr3d-python-latest/fr3d/search/queries/
```

Each `.json` file is one motif query. To run a different set, list the file
names (with or without `.json`) in `query_families`. Examples available in that
folder include sarcin-ricin, kink-turn, GNRA hairpin, internal loop, hairpins,
tandem GU pairs, and junctions.

To run **everything** (the complete FR3D pipeline), just set:

```json
"query_selection": "all"
```

---

## 5. Troubleshooting

- **PyMOL seems to pause on `rmv_db FR3D`:** the search runs synchronously, so
  the window is unresponsive while it works. Only the **first ever run** is
  slow — it builds the persistent annotation cache in
  `output/fr3d_runs/_fr3d_cache/` (including the large reference structures the
  `geometric_*` queries embed). Every run after that reuses the cache and
  finishes in seconds. If a run pauses for minutes, it is building this cache
  for the first time; let it finish once.
- **`all` takes a long time:** it runs every query in the folder; some
  `geometric_*` queries download a reference structure. Use `families` for a
  chosen subset, `all` only for a full analysis.
- **FR3D re-annotates on every run / stays slow:** confirm the
  `output/fr3d_runs/_fr3d_cache/units/` folder is being populated with
  `*_NA.pickle` files after the first run. If it is empty, the interpreter
  could not write there — check folder permissions for `output/fr3d_runs/`.
- **"no usable interpreter" on setup:** run `rmv_setup FR3D` again, or point it
  at a Python that has internet + pip: `rmv_setup FR3D /path/to/python`.
- **A query is skipped:** a name in `query_families` didn't match a file in the
  queries folder, or that query is malformed — the log names it and the rest
  still run.
