# Data Sources

RSMViewer supports **8 data sources** for RNA motif annotations. Each source can be selected
independently or combined for multi-source analysis.

---

## Source Overview

| ID | Name | Type | Coverage | Command |
|----|------|------|----------|---------|
| 1 | RNA 3D Atlas | Local (offline) | 759 PDBs | `rmv_db 1` |
| 2 | Rfam | Local (offline) | 173 PDBs | `rmv_db 2` |
| 3 | BGSU RNA 3D Hub | Online API | ~3000+ PDBs | `rmv_db 3` |
| 4 | Rfam API | Online API | All Rfam families | `rmv_db 4` |
| 5 | FR3D Annotations | User annotations | Custom | `rmv_db 5` |
| 6 | RNAMotifScan (RMS) | User annotations | Custom | `rmv_db 6` |
| 7 | RNAMotifScanX (RMSX) | User annotations | Custom | `rmv_db 7` |
| 8 | NoBIAS | User annotations | Custom | `rmv_db 8` |

---

## Source Configuration

Source IDs are defined in `database/config.py` via `SOURCE_ID_MAP`:

```python
SOURCE_ID_MAP: Dict[int, Dict] = {
    1: {'name': 'RNA 3D Atlas',       'type': 'local', 'subtype': 'atlas',    'command': 'rmv_db 1'},
    2: {'name': 'Rfam',               'type': 'local', 'subtype': 'rfam',     'command': 'rmv_db 2'},
    3: {'name': 'BGSU RNA 3D Hub',    'type': 'web',   'subtype': 'bgsu',     'command': 'rmv_db 3'},
    4: {'name': 'Rfam API',           'type': 'web',   'subtype': 'rfam_api', 'command': 'rmv_db 4'},
    5: {'name': 'FR3D Annotations',   'type': 'user',  'tool': 'fr3d',        'command': 'rmv_db 5'},
    6: {'name': 'RNAMotifScan (RMS)', 'type': 'user',  'tool': 'rms',         'supports_filtering': True, 'command': 'rmv_db 6'},
    7: {'name': 'RNAMotifScanX (RMSX)','type': 'user', 'tool': 'rmsx',        'supports_filtering': True, 'command': 'rmv_db 7'},
    8: {'name': 'NoBIAS',             'type': 'user',  'tool': 'nobias',      'supports_filtering': True, 'command': 'rmv_db 8'},
}
```

Each entry also contains `category`, `mode`, `description`, and `coverage` fields.

---

## Local Sources (1-2)

### Source 1 — RNA 3D Atlas

- **Provider:** `database/atlas_provider.py`
- **Data files:** Bundled JSON in `motif_database/RNA 3D motif atlas/`
- **Motif types:** HL, IL, J3, J4, J5, J6, J7
- **How it works:** Parses pre-packaged JSON files containing loop annotations from the RNA 3D Motif Atlas project. No internet required.

### Source 2 — Rfam

- **Provider:** `database/rfam_provider.py`
- **Data files:** Bundled Stockholm alignment files in `motif_database/Rfam motif database/`
- **Motif types:** GNRA, UNCG, CUYG, K-TURN, T-LOOP, C-LOOP, U-TURN, SARCIN-RICIN, TANDEM-GA, and more (19 types total)
- **How it works:** Parses Stockholm-format seed alignment files to extract motif annotations.

---

## Online Sources (3-4)

### Source 3 — BGSU RNA 3D Hub

- **Provider:** `database/bgsu_api_provider.py`
- **Strategy:** Hybrid HTML scraping + CSV fallback from `rna.bgsu.edu`
- **Coverage:** ~3000+ PDB structures
- **How it works:**
  1. Scrapes `http://rna.bgsu.edu/rna3dhub/nrlist` for loop annotations
  2. Falls back to CSV download if HTML parsing fails
  3. Responses are cached for 30 days via `cache_manager.py`

### Source 4 — Rfam API

- **Provider:** `database/rfam_api_provider.py`
- **Endpoint:** Rfam REST API
- **Coverage:** All Rfam-annotated motif families
- **How it works:** Queries the Rfam REST API for motif family annotations. Cached for 30 days.

---

## User Annotation Sources (5-8)

### Source 5 — FR3D

- **Provider:** `database/user_annotations/user_provider.py`
- **Default directory:** `database/user_annotations/fr3d/`
- **Format:** FR3D output CSV files

### Source 6 — RNAMotifScan (RMS)

- **Provider:** `database/user_annotations/user_provider.py`
- **Default directory:** `database/user_annotations/RNAMotifScan/`
- **Format:** Organized by motif type in subdirectories
- **P-value filtering:** Supported (enabled by default)

### Source 7 — RNAMotifScanX (RMSX)

- **Provider:** `database/user_annotations/user_provider.py`
- **Default directory:** `database/user_annotations/RNAMotifScanX/`
- **Format:** Same structure as RMS with consensus suffixes
- **P-value filtering:** Supported (enabled by default)

### Source 8 — NoBIAS

- **Provider:** `database/user_annotations/user_provider.py`
- **Default directory:** `database/user_annotations/NoBIAS/`
- **P-value filtering:** Supported (enabled by default)

---

## Custom Data Paths

Each user annotation source (5-8) supports a custom data directory. Paths are stored per-source and do not overwrite each other:

```
rmv_db 5 /path/to/fr3d/data       # FR3D custom directory
rmv_db 6 ~/my_rms_data            # RMS custom directory
rmv_db 7 /path/to/rmsx/data       # RMSX custom directory
rmv_db 8 /data/nobias_output      # NoBIAS custom directory
```

Internally, paths are stored in `gui.user_data_paths: Dict[int, str]`.

---

## P-Value Filtering (RMS/RMSX/NoBIAS)

Sources 6, 7, and 8 include P-values in their output. The plugin can filter results based on these values.

```
rmv_db 6 off                     # Disable filtering
rmv_db 6 on                      # Enable filtering
rmv_db 6 SARCIN-RICIN 0.01       # Custom threshold per motif type
rmv_db 7 C-LOOP 0.05             # RMSX custom threshold
rmv_db 8 KINK-TURN 0.02          # NoBIAS custom threshold
```

Custom P-value thresholds are stored in `gui.user_rms_custom_pvalues`, `gui.user_rmsx_custom_pvalues`, and `gui.user_nobias_custom_pvalues` dictionaries.

---

## Source Selection Workflow

```
rmv_fetch 1S72           # Step 1: Load structure
rmv_sources              # Check available sources
rmv_db 3                 # Select BGSU API
rmv_load_motif           # Fetch motif data from BGSU
rmv_summary              # View results
```

Switching sources does not require reloading the PDB structure:

```
rmv_db 1                 # Switch to Atlas
rmv_load_motif           # Fetch from Atlas (PDB stays loaded)
```

---

## Multi-Source Combine

Combine multiple sources in a single command:

```
rmv_db 1 3               # Atlas + BGSU
rmv_db 1 3 4             # Atlas + BGSU + Rfam API
rmv_db 1 3 jaccard_threshold=0.80  # Custom merge threshold
```

See [multi_source_combine.md](multi_source_combine.md) for details on the cascade merge pipeline.

---

## Related Commands

| Command | Purpose |
|---------|---------|
| `rmv_db <N>` | Select data source |
| `rmv_sources` | List all 8 sources |
| `rmv_source info` | Show active source details |
| `rmv_source info <N>` | Show specific source details |
| `rmv_load_motif` | Fetch motif data from selected source |
| `rmv_refresh` | Force API cache refresh |
| `rmv_loaded` | List all loaded PDBs and sources |
