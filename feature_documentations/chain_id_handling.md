# Chain ID Handling

RSMViewer implements a **4-layer chain ID protection system** to ensure consistency between
PDB structure chain IDs in PyMOL and motif annotation data.

---

## Background

PDB/CIF files contain two types of chain identifiers:

| Type | Field | Description | Example |
|------|-------|-------------|---------|
| **auth_asym_id** | Author-assigned | Chain IDs from the depositing authors | `0`, `9`, `A`, `B` |
| **label_asym_id** | System-assigned | Standard chain IDs from the PDB archive | `A`, `AA`, `BA`, `CA` |

By default, RSMViewer uses `auth_asym_id` (the PyMOL default). Motif databases (BGSU, Atlas, Rfam)
also use `auth_asym_id`, so everything matches by default.

---

## The 4-Layer System

### Layer 1 — Startup Lock (`plugin.py`)

At plugin initialization:

```python
cmd.set("cif_use_auth", 1)
```

This locks PyMOL to `auth_asym_id` mode. Without this, PyMOL might use a different
default depending on the user's configuration.

### Layer 2 — Fetch-Time Override (`gui.py`)

When loading with `cif_use_auth=0`:

```
rmv_fetch 1S72 cif_use_auth=0
```

The plugin:
1. Sets `cmd.set("cif_use_auth", 0)` before fetching
2. Stores `gui.cif_use_auth = 0` for later reference
3. All subsequent chain operations use `label_asym_id`

### Layer 3 — CIF File Parsing (`gui.py`)

When `cif_use_auth=0`, the plugin parses the on-disk CIF file's `_atom_site` loop to
build an `auth → label` chain mapping:

```python
# Resulting mapping: {auth_asym_id: label_asym_id}
# Example for 1S72: {'0': 'A', '9': 'B', 'A': 'C', 'B': 'D', ...}
```

This mapping is stored in `gui.auth_to_label_map` and applied during motif data loading
to remap residue chain IDs from auth (used by databases) to label (used by PyMOL).

### Layer 4 — Selection Building (`utils/selectors.py`)

The `SelectionParser` class builds PyMOL selection strings with fallback logic:

```python
# Primary selection: uses "chain" keyword
selection = f"chain {chain_id} and resi {residue_range}"

# If no atoms found, tries "segi" fallback
# (label_asym_id is stored in PyMOL's segi field when cif_use_auth=1)
segi_selection = f"segi {chain_id} and resi {residue_range}"
```

**Verification:** After building a selection, the plugin verifies that atoms were actually
found (`cmd.count_atoms()`). If zero atoms are found with the primary `chain` keyword,
it tries the `segi` keyword as a fallback.

---

## Chain Selection Syntax

PyMOL selections use compact range notation:

```python
# Single chain, continuous range
"chain A and resi 158-164"

# Single chain, multiple ranges
"chain A and resi 158-164+171-178"

# Multi-chain motif
"(chain A and resi 158-164) or (chain B and resi 75-81)"
```

The selection builder in `utils/parser.py` validates chains against the loaded structure:

```python
actual_chains = pymol_cmd.get_chains(structure_name)
if actual_chains and chain not in actual_chains:
    logger.warning(f"Chain '{chain}' not found in '{structure_name}'...")
```

---

## When to Use Label Mode

Use `cif_use_auth=0` when:
- Your user annotations (FR3D/RMS/RMSX/NoBIAS) use `label_asym_id` chain IDs
- You need `label_asym_id` chains for downstream analysis
- Your PDB has unusual auth chain IDs (e.g., numeric: `0`, `9`)

### Example: Auth vs Label

**Default mode** (`cif_use_auth=1`):
```
rmv_fetch 1S72
rmv_chains                   # Chains: 0, 9, A, B, C, D, ...
rmv_summary SARCIN-RICIN     # Chain 0:158-164, 9:75-81
```

**Label mode** (`cif_use_auth=0`):
```
rmv_fetch 1S72 cif_use_auth=0
rmv_chains                   # Chains: A, AA, AB, AC, ... (293 chains)
rmv_summary SARCIN-RICIN     # Chain A:158-164, B:75-81
```

---

## Diagnostics

Use `rmv_chains` to verify chain mode and see loaded chains:

```
rmv_chains
```

Output:
```
Structure: 1S72  |  cif_use_auth = 0 (label_asym_id)  |  Chains: 293
Label chains:  A AA AB AC AD AE ...
```

---

## Related Commands

| Command | Purpose |
|---------|---------|
| `rmv_fetch <PDB> cif_use_auth=0` | Load with label chain IDs |
| `rmv_chains` | Show chain ID mode and loaded chains |
| `rmv_loaded` | List loaded PDBs with chain mode info |
