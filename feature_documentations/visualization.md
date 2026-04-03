# Visualization

RSMViewer provides two primary visualization commands: `rmv_show` (creates PyMOL objects) and `rmv_view` (lightweight in-place highlighting).

---

## rmv_show — Full Visualization

Creates separate PyMOL objects for motif instances, applies coloring, and optionally zooms to specific instances.

### Usage

```
rmv_show HL                # Show all hairpin loops (creates HL_ALL_S3 object)
rmv_show GNRA 1            # Zoom to GNRA instance #1 (creates GNRA_1_S3)
rmv_show HL 1,2,3          # Show specific instances
rmv_show ALL               # Show all loaded motif types
rmv_show K-TURN, padding=3 # Expand view ±3 residues for context
```

### How It Works

1. Resolves the motif type name (handles multi-word names like `4-WAY JUNCTION (J4)`)
2. Looks up motif instances from loaded data
3. Creates a PyMOL selection using chain + residue range notation
4. Creates a new PyMOL object with a source suffix (e.g., `HL_ALL_S3` for source 3)
5. Applies motif-specific color from `colors.py`
6. Sets non-motif residues to background color (default: `gray80`)
7. If an instance number is given: zooms camera and prints residue details

### Source-Tagged Objects

When working with multiple sources, each source gets its own tagged objects:

```
rmv_db 3
rmv_show HL            # Creates: HL_ALL_S3, HL_1_S3, HL_2_S3, ...
rmv_db 1
rmv_show HL            # Creates: HL_ALL_S1, HL_1_S1, ... (S3 objects persist)
```

This allows direct source comparison using native PyMOL commands:

```
align HL_3_S3, HL_3_S1
```

### Padding

The `padding=N` parameter expands the visualized residue range by ±N positions on each side. This is useful for seeing structural context around a motif.

```
rmv_show GNRA 1, padding=5     # Instance 1 with 5 flanking residues
rmv_show K-TURN, padding=3     # All K-turns with 3 flanking residues
```

Padding only affects **visualization** (coloring and PyMOL objects). Summary tables, instance details, and stored data are unchanged.

---

## rmv_view — Lightweight Highlighting

Highlights and zooms to motif regions directly on the base structure **without creating new PyMOL objects**. Faster than `rmv_show` for quick exploration.

### Usage

```
rmv_view HL              # Highlight all hairpin loop residues
rmv_view HL 1            # Highlight and zoom to HL instance 1
rmv_view GNRA 3          # Highlight and zoom to GNRA instance 3
rmv_view MOTIF           # Auto-color all motif regions on structure
```

### How It Works

The command is implemented in `gui.py` (`view_motif` function):

1. If called with `MOTIF` as the type → calls `_auto_color_motifs_on_structure()` to highlight all motif regions
2. Otherwise, resolves the motif type and optional instance numbers
3. With instance numbers → calls `viz_manager.view_motif_instance()` per instance (highlight + zoom)
4. Without instance numbers → calls `viz_manager.view_motif_type()` for all instances

### When to Use rmv_view vs rmv_show

| Feature | `rmv_view` | `rmv_show` |
|---------|-----------|-----------|
| Creates PyMOL objects | No | Yes |
| Speed | Faster | Slower (object creation) |
| Supports `rmv_save` | No (no objects) | Yes |
| Supports `rmv_super`/`rmv_align` | No | Yes |
| Good for exploration | Yes | Yes |
| Good for publication figures | No | Yes |

---

## rmv_toggle — Visibility Control

Toggles visibility of motif type objects in the PyMOL viewport.

```
rmv_toggle HL off        # Hide hairpin loop objects
rmv_toggle HL on         # Show them again
```

---

## rmv_bg_color — Background Color

Changes the color of non-motif residues (the structural backbone).

```
rmv_bg_color white       # White backbone
rmv_bg_color lightgray   # Light gray
rmv_bg_color gray80      # Default
```

---

## Color System

### Default Colors (`colors.py`)

Colors are defined as RGB tuples normalized to 0–1:

```python
MOTIF_COLORS = {
    'HL': (1.0, 0.0, 0.0),      # Red
    'IL': (0.0, 1.0, 1.0),      # Cyan
    'J3': (1.0, 1.0, 0.0),      # Yellow
    'J4': (1.0, 0.0, 1.0),      # Magenta
    'J5': (0.0, 1.0, 0.0),      # Green
    'J6': (1.0, 0.5, 0.0),      # Orange
    'J7': (0.0, 0.0, 1.0),      # Blue
    'GNRA': (0.0, 0.8, 0.4),    # Teal green
    'K-TURN': (0.2, 0.6, 1.0),  # Bright blue
    'SARCIN-RICIN': (0.8, 0.0, 0.2),  # Dark red
    'T-LOOP': (1.0, 0.4, 0.7),  # Pink
    ...
}
```

Unknown motif types get `DEFAULT_COLOR = (1.0, 0.5, 0.0)` (bright orange).

### Custom Colors

Override any motif color at runtime:

```
rmv_color HL blue        # Change HL to blue
rmv_color GNRA red       # Change GNRA to red
rmv_colors               # View current color assignments
```

Custom colors are stored in `CUSTOM_COLORS` and take precedence over `MOTIF_COLORS`.

---

## Typical Visualization Workflow

```
rmv_fetch 1S72
rmv_db 3
rmv_load_motif

rmv_summary              # See what motif types are available
rmv_view GNRA 1          # Quick peek at GNRA instance 1
rmv_show GNRA            # Create persistent objects for all GNRA
rmv_show GNRA 1          # Zoom to instance 1
rmv_color GNRA red       # Change GNRA color
rmv_bg_color white       # White backbone
rmv_save GNRA sticks     # Save images
```
