"""
RSMViewer - Logging Module
Provides unified logging for the plugin with optional PyMOL console output.
"""

import sys
from datetime import datetime


class PluginLogger:
    """Simple logging system for the RSMViewer plugin."""
    
    def __init__(self, use_pymol_console=False):
        """
        Initialize logger.
        
        Args:
            use_pymol_console (bool): If True, also print to PyMOL console
        """
        self.use_pymol_console = use_pymol_console
        self.log_file = None
    
    def set_log_file(self, filepath):
        """Set optional log file path."""
        self.log_file = filepath
    
    def _format_message(self, level, message):
        """Format log message with timestamp."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        return f"[{timestamp}] [{level}] {message}"
    
    def _write(self, formatted_msg):
        """Write message to appropriate outputs."""
        print(formatted_msg)
        
        if self.log_file:
            try:
                with open(self.log_file, 'a') as f:
                    f.write(formatted_msg + '\n')
            except Exception as e:
                print(f"Warning: Could not write to log file: {e}")
    
    def info(self, message):
        """Log info level message."""
        self._write(self._format_message("INFO", message))
    
    def warning(self, message):
        """Log warning level message."""
        self._write(self._format_message("WARNING", message))
    
    def error(self, message):
        """Log error level message."""
        self._write(self._format_message("ERROR", message))
    
    def debug(self, message):
        """Log debug level message."""
        self._write(self._format_message("DEBUG", message))
    
    def success(self, message):
        """Log success level message."""
        self._write(self._format_message("SUCCESS", message))


# Global logger instance
_logger = None


def get_logger():
    """Get or create global logger instance."""
    global _logger
    if _logger is None:
        _logger = PluginLogger(use_pymol_console=True)
    return _logger


def initialize_logger(use_pymol_console=False, log_file=None):
    """Initialize global logger."""
    global _logger
    _logger = PluginLogger(use_pymol_console=use_pymol_console)
    if log_file:
        _logger.set_log_file(log_file)
    return _logger
