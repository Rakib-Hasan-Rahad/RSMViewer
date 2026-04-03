# Image Export & Structure Export

RSMViewer supports exporting motif instances as **PNG images** (via `image_saver.py`) and
as standalone **mmCIF structure files** (via `structure_exporter.py`).

---

## PNG Image Export

### Commands

```
rmv_save ALL                   # All motifs, default representation (cartoon)
rmv_save ALL sticks            # All motifs, sticks representation
rmv_save HL                    # All hairpin loops
rmv_save HL 3 spheres          # HL instance 3 as spheres
rmv_save current               # Current PyMOL viewport (high-res)
rmv_save current my_view       # Custom filename
```

### Supported Representations

The `MotifImageSaver` class in `image_saver.py` supports these representations:

| Representation | Description |
|---------------|-------------|
| `cartoon` | Default — shows secondary structure |
| `sticks` | Ball-and-stick atomic detail |
| `spheres` | Space-filling model |
| `ribbon` | Simplified ribbon trace |
| `lines` | Wire frame |
| `licorice` | Thick stick model |
| `surface` | Molecular surface |
| `cartoon+sticks` | Combined view (split on `+`) |

### How It Works

1. **Isolate:** Creates a temporary PyMOL object containing only the motif residues
2. **Style:** Applies the selected representation and motif-type coloring
3. **Frame:** Zooms camera to the isolated motif
4. **Render:** Captures at 800×600 resolution
5. **Clean:** Removes the temporary object

For `rmv_save current`:
- Captures the current viewport as-is (no modifications)
- Resolution: 2400×1800 at 300 DPI (high-res)

### Output Folder

```
motif_images/
└── <pdb_id>/
    ├── <MOTIF_TYPE>/
    │   ├── <TYPE>-<NO>-<chain>_<residues>.png
    │   └── ...
    └── current_view_<timestamp>.png
```

### Key Functions

```python
class MotifImageSaver:
    def save_instance_image(self, motif_folder, instance_no, motif_type,
                           motif_details, structure_name, representation='cartoon') -> bool

    def save_all_motifs(self, loaded_motifs, structure_name, pdb_id,
                       output_base_dir=None, representation='cartoon') -> Dict

    def save_motif_type_images(self, loaded_motifs, motif_type, structure_name,
                              pdb_id, output_base_dir=None, representation='cartoon') -> Dict
```

---

## mmCIF Structure Export

### Commands

```
rmv_save ALL cif               # Export ALL motif instances as mmCIF
rmv_save HL cif                # Export all hairpin loops as mmCIF
rmv_save HL 3 cif              # Export HL instance #3 as mmCIF
```

### How It Works

The `MotifStructureExporter` class in `structure_exporter.py` extracts motif instances
as standalone `.cif` files using **original coordinates from the on-disk CIF file**
(not PyMOL's internal coordinates, which may be slightly modified during loading).

#### Algorithm

1. **Find CIF:** `_find_cif_file(pdb_id)` locates the original `.cif` on disk via PyMOL's `fetch_path`
2. **Build residue set:** Creates `{(chain_str, resi_str)}` tuples from motif instance details
3. **Parse headers:** `_parse_atom_site_header()` finds `auth_asym_id` (chain) and `auth_seq_id` (residue number) column indices in the `_atom_site` loop
4. **Filter atoms:** Iterates data rows, tokenizes each, checks if `(tokens[auth_asym_col], tokens[auth_seq_col])` is in the residue set
5. **Write CIF:** Outputs a minimal mmCIF file containing:
   - Metadata blocks (`_cell`, `_symmetry`, `_entity`, etc.) copied from source
   - `loop_` + `_atom_site` column headers + matched atom rows only

#### Why Not Use PyMOL's `cmd.save()`?

PyMOL may modify atomic coordinates slightly during loading (symmetry expansion, origin
shifting). The structure exporter bypasses PyMOL entirely and reads/writes the raw CIF text
to preserve original deposited coordinates.

### Output Folder

```
motif_structures/
└── <pdb_id>/
    ├── <MOTIF_TYPE>/
    │   ├── <TYPE>-1-<chain>_<residues>.cif
    │   ├── <TYPE>-2-<chain>_<residues>.cif
    │   └── ...
    └── ...
```

Each `.cif` file is a self-contained mmCIF that can be loaded independently:

```
cmd.load("motif_structures/1s72/HL/HL-3-0_55-64.cif")
```

---

## Typical Export Workflow

```
rmv_fetch 1S72
rmv_db 3
rmv_load_motif

rmv_show GNRA                 # Visualize first
rmv_save GNRA sticks          # Save PNG images
rmv_save GNRA cif             # Export mmCIF structures
rmv_save current              # Save current viewport
```

---

## Related Commands

| Command | Purpose |
|---------|---------|
| `rmv_show <TYPE>` | Visualize before saving |
| `rmv_summary <TYPE>` | Check instance count before batch export |
| `rmv_super <TYPE>` | Superimpose before saving superimposition view |
