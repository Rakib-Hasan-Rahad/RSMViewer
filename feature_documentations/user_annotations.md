# User Annotations

RSMViewer supports loading motif annotations from four external tools: **FR3D** (source 5),
**RNAMotifScan** (source 6), **RNAMotifScanX** (source 7), and **NoBIAS** (source 8).

---

## Supported Tools

| ID | Tool | Provider | P-Value Filtering |
|----|------|----------|-------------------|
| 5 | FR3D | `database/user_annotations/user_provider.py` | No |
| 6 | RNAMotifScan (RMS) | Same | Yes |
| 7 | RNAMotifScanX (RMSX) | Same | Yes |
| 8 | NoBIAS | Same | Yes |

---

## Directory Setup

Place annotation files in the default directories:

```
rsmviewer/database/user_annotations/
├── fr3d/                        # FR3D output files
│   └── <PDB>_motifs.csv
├── RNAMotifScan/                # RMS files (by motif type)
│   ├── c_loop/
│   │   └── Res_<pdb_id>
│   ├── sarcin_ricin/
│   │   └── Res_<pdb_id>
│   ├── Kturn/
│   │   └── Res_<pdb_id>
│   └── ...
├── RNAMotifScanX/               # RMSX files (same structure)
│   ├── c-loop_consensus/
│   │   └── Res_<pdb_id>
│   ├── k-turn_consensus/
│   │   └── Res_<pdb_id>
│   └── ...
└── NoBIAS/                      # NoBIAS output files
    └── <pdb_id>_<motif>_nobias.txt
```

---

## Loading User Annotations

```
rmv_fetch 1S72

rmv_db 5               # FR3D
rmv_load_motif

rmv_db 6               # RNAMotifScan (RMS)
rmv_load_motif

rmv_db 7               # RNAMotifScanX (RMSX)
rmv_load_motif

rmv_db 8               # NoBIAS
rmv_load_motif
```

---

## Custom Data Paths (Per-Source)

Each source stores its own custom path independently:

```
rmv_db 5 /path/to/fr3d/data       # FR3D custom directory
rmv_db 6 ~/my_rms_data            # RMS custom directory
rmv_db 7 /path/to/rmsx/data       # RMSX custom directory
rmv_db 8 /data/nobias_output      # NoBIAS custom directory
```

Paths are stored in `gui.user_data_paths: Dict[int, str]`. Setting a path for one source does **not** affect other sources.

---

## P-Value Filtering (RMS/RMSX/NoBIAS)

Sources 6, 7, and 8 include P-values in their output. The plugin filters results based on these values.

### Enable/Disable Filtering

```
rmv_db 6 off                     # Disable RMS filtering (show ALL results)
rmv_db 6 on                      # Re-enable RMS filtering

rmv_db 7 off                     # Disable RMSX filtering
rmv_db 8 off                     # Disable NoBIAS filtering
```

### Custom Thresholds Per Motif Type

```
rmv_db 6 SARCIN-RICIN 0.01       # RMS: only SARCIN-RICIN with P ≤ 0.01
rmv_db 7 C-LOOP_CONSENSUS 0.05   # RMSX custom threshold
rmv_db 8 KINK-TURN 0.02          # NoBIAS custom threshold
```

### How Filtering Works

When filtering is enabled:
- Each motif instance's P-value is compared against the threshold
- Only instances with P-value ≤ threshold are included
- Default thresholds are built into the plugin
- Custom thresholds override the defaults for specific motif types

Filtering state is stored in:
- `gui.user_rms_filtering_enabled` (bool)
- `gui.user_rmsx_filtering_enabled` (bool)
- `gui.user_nobias_filtering_enabled` (bool)

Custom thresholds are stored in:
- `gui.user_rms_custom_pvalues` (dict)
- `gui.user_rmsx_custom_pvalues` (dict)
- `gui.user_nobias_custom_pvalues` (dict)

---

## Format Parsers

The `database/user_annotations/converters.py` module parses each tool's output format:

- **FR3D:** CSV format with chain, residue number, and nucleotide columns
- **RMS/RMSX:** Tab/space-delimited files organized by motif type in subdirectories
- **NoBIAS:** Custom text format with P-value scores

The `user_provider.py` module loads and converts these formats into the unified `MotifInstance` / `MotifType` data structures used throughout the plugin.

---

## Combining User Annotations with Other Sources

User annotation sources can be combined with local/API sources:

```
rmv_db 1 6               # Atlas + RMS
rmv_db 3 7               # BGSU + RMSX
rmv_db 1 3 8             # Atlas + BGSU + NoBIAS
```

The cascade merge pipeline handles deduplication across all sources. See [multi_source_combine.md](multi_source_combine.md) for details.

---

## Legacy Command: rmv_user

The older `rmv_user` command is still available:

```
rmv_user fr3d 1S72
rmv_user rnamotifscan 1A00
rmv_user list                # Show available annotation files
```

The preferred approach is `rmv_db <N>` + `rmv_load_motif`.

---

## Related Commands

| Command | Purpose |
|---------|---------|
| `rmv_db <N>` | Select user annotation source (5-8) |
| `rmv_db <N> off/on` | Toggle P-value filtering |
| `rmv_db <N> <MOTIF> <P>` | Set custom P-value threshold |
| `rmv_db <N> /path` | Set custom data directory |
| `rmv_load_motif` | Load annotations from selected source |
| `rmv_summary` | View loaded annotation results |
| `rmv_show <TYPE>` | Visualize annotated motifs |
| `rmv_view <TYPE>` | Quick highlight annotated motifs |
