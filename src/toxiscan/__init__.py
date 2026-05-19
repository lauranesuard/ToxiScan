"""
ToxiScan - A Python package to evaluate molecular toxicity.

This package provides tools to:
- Convert molecule names to SMILES via PubChem
- Detect toxic functional groups (toxicophores)
- Calculate a toxicity score
- Visualize molecules with highlighted toxic fragments
"""

from .molecule import get_smiles, get_molecule_info
from .toxicophores import find_toxicophores
from .scoring import (
    remove_redundant_toxicophores,
    count_toxicophores,
    toxicity_approximation,
    compute_properties,
    interpret_properties,
    properties_toxicity,
)
from .visualization import draw_molecule, draw_molecule_3d

__all__ = [
    "get_smiles",
    "get_molecule_info",
    "find_toxicophores",
    "remove_redundant_toxicophores",
    "count_toxicophores",
    "toxicity_approximation",
    "compute_properties",
    "interpret_properties",
    "properties_toxicity",
    "draw_molecule",
    "draw_molecule_3d",
]