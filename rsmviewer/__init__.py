"""RSMViewer: a four-source, motif-ID PyMOL plugin."""

from pathlib import Path as _Path


# PyMOL's startup loader may execute this directory's __init__.py as a module
# without assigning package search paths; provide one for relative imports.
if "__path__" not in globals():
    __path__ = [str(_Path(__file__).resolve().parent)]

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
