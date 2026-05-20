import pytest
import os
from src.toxiscan.visualization import draw_molecule, draw_molecule_3d

def test_draw_molecule_invalid_smiles():
    with pytest.raises(ValueError):
        draw_molecule("INVALIDE", {})

def test_draw_molecule():
    # aspirine sans toxicophores
    draw_molecule("CC(=O)Oc1ccccc1C(=O)O", {})
    assert os.path.exists("molecule.svg")

def test_draw_molecule_3d_invalid_smiles():
    with pytest.raises(ValueError):
        draw_molecule_3d("INVALIDE", {})

def test_draw_molecule_3d():
    # formaldéhyde avec un aldéhyde détecté
    result = draw_molecule_3d("C=O", {"Aldehyde": [0, 1]})
    assert isinstance(result, str)  
    assert "<div" in result        

