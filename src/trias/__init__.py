"""TRiAs: Tandem Repeat analysis toolkit.

The version lives here and nowhere else; pyproject.toml reads it via
[tool.hatch.version]. Keep this module import-free so that `from . import
__version__` in cli.py cannot create a circular import.
"""

__version__ = "0.2.0"

__all__ = ["__version__"]
