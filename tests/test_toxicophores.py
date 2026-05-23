from toxiscan.toxicophores import find_toxicophores 
import pytest

def test_find_toxicophores_zero():
    """
    Test that find_toxicophores returns an empty dictionary
    for aspirin, which contains no known toxicophores.
    """
    result = find_toxicophores("CC(=O)Oc1ccccc1C(=O)O")
    assert result == {}

def test_find_toxicophores_non_zero():
    """
    Test that find_toxicophores detects the aldehyde group
    in formaldehyde.
    """
    result = find_toxicophores("C=O")
    assert "Aldehyde" in result 

def test_find_toxicophores_error():
    """
    Test that find_toxicophores raises a ValueError
    when given an invalid SMILES string.
    """
    with pytest.raises(ValueError):
        find_toxicophores("INVALIDES")

def test_find_toxicophores_nitro():
    """
    Test that find_toxicophores correctly detects the nitro group
    in nitrobenzene.
    """
    result = find_toxicophores("c1ccc([N+](=O)[O-])cc1")
    assert "Nitro group" in result


def test_find_toxicophores_returns_atom_indices():
    """
    Test that find_toxicophores returns atom indices as lists
    for each detected toxicophore.
    """
    result = find_toxicophores("c1ccc([N+](=O)[O-])cc1")
    for name, indices in result.items():
        assert isinstance(indices, list)
        assert len(indices) > 0


