# User Annotations Feature

This directory contains support for loading motif annotations from external tools like FR3D and RNAMotifScan.

## Directory Structure

```
user_annotations/
├── fr3d/              # FR3D output files
├── rnamotifscan/      # RNAMotifScan output files
├── converters.py      # Format converters
└── user_provider.py   # User annotation provider
```

## Supported Tools

### FR3D (Fidelity of RNA 3D Structure)
- **Format**: CSV (comma-separated)
- **Columns**: Motif order, Motif type, Resolution, Positions, Sequence, cWW, Description
- **Position Format**: `PDB_ID|Chain|Model|Start-End` (e.g., `1S72|1|0|13-530`)
- **Expected File Names**: `{pdb_id}_motifs.csv` (e.g., `1s72_motifs.csv`)

#### Example FR3D Usage:
1. Generate motif annotation using FR3D tool (external)
2. Save as CSV file: `1s72_motifs.csv`
3. Place in: `database/user_annotations/fr3d/`
4. In PyMOL: `rmv_user fr3d 1S72`
5. Then use: `rmv_summary`, `rmv_show`, etc.

### RNAMotifScan
- **Format**: CSV or TSV
- **Columns**: Motif_Name, Start, End, Sequence, Chain, Score, etc.
- **Expected File Names**: `{pdb_id}.csv` or `{pdb_id}.tsv`

#### Example RNAMotifScan Usage:
1. Generate motif annotation using RNAMotifScan tool (external)
2. Save as CSV/TSV file
3. Place in: `database/user_annotations/rnamotifscan/`
4. In PyMOL: `rmv_user rnamotifscan 1A00`
5. Then use standard commands

## PyMOL Commands

### Load User Annotations
```pymol
# Show help
rmv_user

# Load FR3D annotations for a structure
rmv_user fr3d 1S72

# Load RNAMotifScan annotations
rmv_user rnamotifscan 1A00

# List all available annotation files
rmv_user list
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

User annotation tools (FR3D, RNAMotifScan) are external tools.
Please cite appropriately:
- FR3D: https://www.bgsu.edu/research/rna/software/fr3d.html
- RNAMotifScan: http://bioinformatics.bc.edu/rnamotif/
