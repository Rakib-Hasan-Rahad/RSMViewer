"""RSMViewer — RNA Structural Motif Viewer

PyMOL plugin for visualizing RNA structural motifs from multiple databases.

This package provides:
- 8 data sources (Atlas, Rfam, BGSU API, Rfam API, FR3D, RMS, RMSX, NoBIAS)
- Multi-source combine with cascade merge
- Label chain ID support (cif_use_auth=0)
- Custom data path support for user annotation sources (FR3D, RMS, RMSX, NoBIAS)
- High-resolution image export with 8 representation options
- Base-pair visualization with Leontis-Westhof labels
- Medoid-based structural superimposition
- 24 PyMOL commands (rmv_help for full reference)

Author: CBB LAB @Rakib Hasan Rahad
Version: 1.0.0
Last Updated: 02 April 2026
"""

# NOTE:
# Keep this package importable outside PyMOL.
# PyMOL injects/provides the `pymol` module at runtime; importing it in a
# regular Python interpreter (e.g., during CLI tests) will fail.

__version__ = '1.0.0'
__author__ = 'CBB LAB KU @Rakib Hasan Rahad'


def __getattr__(name):
    """Lazy imports so non-PyMOL environments can import this package."""
    if name == '__init_plugin__':
        from .plugin import __init_plugin__ as value
        return value
    if name in {'get_gui', 'initialize_gui'}:
        from . import gui as _gui
        return getattr(_gui, name)
    if name == 'VisualizationManager':
        from .loader import VisualizationManager as value
        return value
    raise AttributeError(name)


__all__ = ['__init_plugin__', 'get_gui', 'initialize_gui', 'VisualizationManager']
