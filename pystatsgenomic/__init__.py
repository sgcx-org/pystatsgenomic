"""
PyStatsGenomic: Genomics and computational-biology statistical computing.

Part of the PyStatistics open-core ecosystem (alongside ``pystatistics`` and
``pystatsbio``). Provides genomic / statistical-genetics methods built on the
general statistical layer.

Usage:
    from pystatsgenomic import hwe
    result = hwe.hardy_weinberg(298, 489, 213)
"""

__version__ = "0.1.0"
__author__ = "Hai-Shuo"
__email__ = "contact@sgcx.org"

from pystatsgenomic import hwe

__all__ = [
    "__version__",
    "hwe",
]
