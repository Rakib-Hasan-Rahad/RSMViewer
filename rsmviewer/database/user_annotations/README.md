# User Annotations Feature

This directory contains support for loading motif annotations from FR3D, RNAMotifScan (RMS), RNAMotifScanX (RMSX), and NoBIAS.

## Directory Structure

```
user_annotations/
├── fr3d/              # FR3D output files
├── RNAMotifScan/      # RNAMotifScan output files
├── RNAMotifScanX/     # RNAMotifScanX output files
├── NoBIAS/            # NoBIAS output files
├── converters.py      # Format converters
└── user_provider.py   # User annotation provider
```

## Supported Tools

### File Encoding Requirement

All user annotation files must be UTF-8 encoded.

- FR3D CSV files
- RNAMotifScan result files
- RNAMotifScanX result files
- NoBIAS text/log files

If loading fails with decode errors, re-save the file as UTF-8 and try again.

### FR3D (Fidelity of RNA 3D Structure)
- **Format**: CSV (motif table) or FR3D pairwise TXT output
- **Columns**: Motif order, Motif type, Resolution, Positions, Sequence, cWW, Description
- **Position Format**: `PDB_ID|Chain|Model|Start-End` (e.g., `1S72|1|0|13-530`)
- **Expected File Names**: `{pdb_id}_motifs.csv` (e.g., `1s72_motifs.csv`)

FR3D pairwise TXT files should look like:

```text
1S72|1|A|G|71	cWW	1S72|1|A|C|83	0
```

and can be named like `1S72_basepair.txt`.

#### Example FR3D Usage:
1. Select Source 5 in PyMOL: `rmv_db 5`
2. Load a structure: `rmv_fetch 1S72`
3. Load motifs: `rmv_load_motif`
4. Explore: `rmv_summary`, `rmv_show`, etc.

#### FR3D wrapper commands in RSMViewer:

```pymol
rmv_fr3d status
rmv_fr3d doctor
rmv_fr3d setup
rmv_fr3d refresh [PDB_ID]
rmv_fr3d config /absolute/path/to/fr3d-python
rmv_fr3d run 1S72
rmv_fr3d run_current
```

`rmv_fr3d run` executes the local Source-5 FR3D pipeline and loads source 5 motifs automatically.

### RNAMotifScan
- **Format**: CSV or TSV
- **Columns**: Motif_Name, Start, End, Sequence, Chain, Score, etc.
- **Expected File Names**: `{pdb_id}.csv` or `{pdb_id}.tsv`

#### Example RNAMotifScan Usage:
1. Place files in `database/user_annotations/RNAMotifScan/`
2. In PyMOL: `rmv_db 6`
3. Run: `rmv_fetch 1A00`
4. Load motifs: `rmv_load_motif`
5. Then use standard commands

### RNAMotifScanX (RMSX)

- **Directory**: `database/user_annotations/RNAMotifScanX/`
- **Typical file**: `result_0_100_withbs.log` inside `*_consensus/` folders
- **Preferred load path**: `rmv_db 7`, then `rmv_fetch`, then `rmv_load_motif`

### NoBIAS

- **Directory**: `database/user_annotations/NoBIAS/`
- **Preferred load path**: `rmv_db 8`, then `rmv_fetch`, then `rmv_load_motif`

## PyMOL Commands

### Load User Annotations
```pymol
# Show help
rmv_user

# Load FR3D annotations for a structure
rmv_user fr3d 1S72

# Load RNAMotifScan annotations
rmv_user rnamotifscan 1A00

# Load RNAMotifScanX annotations
rmv_user rnamotifscanx 1A00

# List all available annotation files
rmv_user list
```

Preferred runtime workflow:

```pymol
rmv_db 5
rmv_fetch 1S72
rmv_load_motif
```

### Work with User Annotations
Once loaded, use standard commands:
```pymol
rmv_summary              # Show all motifs from annotation file
rmv_summary HL           # Show hairpin loop instances
rmv_show HL              # Render hairpin loops on structure
rmv_show HL 1            # View specific instance
rmv_colors               # Show color legend
```

## File Format Reference

### FR3D CSV Format
```csv
Motif order,Motif type,Resolution,Positions,Sequence,cWW,Description
1,Hairpin,NA,"1S72|1|0|13-530","GCCAGCUGGUUGC...",278,"Hairpin with 10 base pairs"
2,Hairpin,NA,"1S72|1|0|27-516","UGCGGCUCAGGGC...",258,"Hairpin with 10 base pairs"
```

Key fields:
- **Motif type**: Used as motif category (e.g., "Hairpin", "Internal loop", "Bulge")
- **Positions**: Format is critical - must be `PDB_ID|Chain|Model|Start-End`
- **Description**: Optional metadata about the motif

### RNAMotifScan Format
```csv
Motif_Name,Start,End,Sequence,Chain,Score
HL,10,50,CGAA...,A,0.95
IL,100,150,AUAA...,A,0.87
```

Key fields:
- **Motif_Name**: Used as motif category (e.g., "HL", "IL", "GNRA")
- **Start/End**: Residue position range
- **Chain**: Chain identifier
- **Score**: Optional quality metric

## Implementation Details

### Converter Architecture
Each converter follows this pattern:

```python
# 1. Parse tool-specific format
# 2. Extract motif name, residues, and metadata
# 3. Convert to standard MotifInstance format
# 4. Return dict: {motif_type: [MotifInstance, ...]}
```

### Data Flow
```
External Tool Output
    ↓
Tool-Specific Converter (FR3DConverter, RNAMotifScanConverter)
    ↓
MotifInstanceSimple objects
    ↓
UserAnnotationProvider converts to standard format
    ↓
Standard MotifInstance objects
    ↓
GUI can display like any other data source
```

## Adding Support for New Tools

To add support for a new tool:

1. Create a new converter class in `converters.py`:
   ```python
   class NewToolConverter:
       @staticmethod
       def convert_file(csv_path: str) -> Dict[str, List[MotifInstanceSimple]]:
           # Parse format
           # Return {motif_type: [instances]}
   ```

2. Update `UserAnnotationProvider._load_file()` to use it:
   ```python
   if tool_name.lower() == 'newtool':
       return NewToolConverter.convert_file(str(file_path))
   ```

3. Create folder: `user_annotations/newtool/`

4. Add to help text and documentation

## Troubleshooting

### "No annotation files found"
- Check file is in correct directory
- Verify file naming convention (e.g., `1s72_motifs.csv` for PDB ID 1S72)
- Run `rmv_user list` to see available files

### Motifs not displaying correctly
- Verify Position format for FR3D (must be `PDB_ID|Chain|Model|Start-End`)
- Check that residue numbers match your PDB file
- Ensure chain IDs are correct

### Import errors
- Ensure CSV file is UTF-8 encoded
- Check for special characters in file
- Verify column headers match expected format

## Example Workflow

```pymol
# Setup: Load PDB structure
pymol> fetch 1S72

# Load FR3D annotations (pre-generated externally)
pymol> rmv_user fr3d 1S72

# Explore motifs in console
pymol> rmv_summary
# Shows: Available motifs: Hairpin (45 instances), ...

# View specific motif type
pymol> rmv_summary Hairpin
# Shows: Table of all 45 hairpin instances

# Render on structure
pymol> rmv_show Hairpin

# View specific instance
pymol> rmv_show Hairpin 1
# Shows: Details of hairpin instance #1
```

## Performance Notes

- User annotations are loaded from local files (fast)
- No API calls required
- Data is cached in memory during session
- Large annotation files (~1000+ motifs) work fine

## License & Attribution

User annotation tools (FR3D, RNAMotifScan, RNAMotifScanX, NoBIAS) are external tools.
Please cite appropriately:
- FR3D: https://www.bgsu.edu/research/rna/software/fr3d.html
- RNAMotifScan: http://bioinformatics.bc.edu/rnamotif/

For full FR3D runtime architecture, setup, and troubleshooting, see [fr3d_runtime/README.md](../../tools/fr3d_runtime/README.md). For RMSX runtime details, see [rmsx_runtime/README.md](../../tools/rmsx_runtime/README.md).
