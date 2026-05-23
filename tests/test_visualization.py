import pytest
import os
from toxiscan.visualization import draw_molecule, draw_molecule_3d

def test_draw_molecule_invalid_smiles():
    """
    Test that draw_molecule raises a ValueError
    when given an invalid SMILES string.
    """
    with pytest.raises(ValueError):
        draw_molecule("INVALIDE", {})

def test_draw_molecule():
    """
    Test that draw_molecule creates a molecule.svg file
    for a valid SMILES string without toxicophores.
    """
    # aspirine sans toxicophores
    draw_molecule("CC(=O)Oc1ccccc1C(=O)O", {})
    assert os.path.exists("molecule.svg")

def test_draw_molecule_3d_invalid_smiles():
    """
    Test that draw_molecule_3d raises a ValueError
    when given an invalid SMILES string.
    """
    with pytest.raises(ValueError):
        draw_molecule_3d("INVALIDE", {})

def test_draw_molecule_3d():
    """
    Test that draw_molecule_3d returns a valid HTML string
    containing a 3D viewer for formaldehyde with highlighted aldehyde group.
    """
    # formaldéhyde avec un aldéhyde détecté
    result = draw_molecule_3d("C=O", {"Aldehyde": [0, 1]})
    assert isinstance(result, str)  
    assert "<div" in result        

