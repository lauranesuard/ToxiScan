from src.toxiscan.toxicophores import find_toxicophores 
import pytest

def test_find_toxicophores_zero():
    result = find_toxicophores("CC(=O)Oc1ccccc1C(=O)O")
    assert result == {}

def test_find_toxicophores_non_zero():
    result = find_toxicophores("C=O")
    assert "Aldehyde" in result 

def test_find_toxicophores_error():
    with pytest.raises(ValueError):
        find_toxicophores("INVALIDES")



