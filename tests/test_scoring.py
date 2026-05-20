from src.toxiscan.scoring import remove_redundant_toxicophores, count_toxicophores, compute_properties, interpret_properties, properties_toxicity, toxicity_approximation
import pytest

def test_remove_redundant_no_redundancy():
    input_dict = {
        "Epoxide": (1, 2, 3),
        "Nitro group": (4, 5, 6),
    }
    result = remove_redundant_toxicophores(input_dict)
    assert "Epoxide" in result
    assert "Nitro group" in result

def test_remove_redundant_with_redundancy():
    input_dict = {
        "Epoxide": (1, 2),
        "Nitro group": (1, 2, 3, 4),
    }
    result = remove_redundant_toxicophores(input_dict)
    assert "Nitro group" in result
    assert "Epoxide" not in result

def test_count_toxicophores_zero():
    result = count_toxicophores("CC(=O)Oc1ccccc1C(=O)O")  # aspirine
    assert result == 0

def test_count_toxicophores_nonzero():
    result = count_toxicophores("C=O")  # formaldéhyde, 1 aldéhyde
    assert result == 1

def test_toxicity_approximation_zero():
    result = toxicity_approximation("CC(=O)Oc1ccccc1C(=O)O")  
    assert result == 0.0

# on re test la fonction dans le cas où un toxicophore est détecté
def test_toxicity_approximation_nonzero():
    result = toxicity_approximation("C=O")  # formaldéhyde, aldéhyde
    assert result == 1.0

def test_compute_properties_fail():
    with pytest.raises(ValueError): # pour signifier je m'attends à une ValueError -> normal, le code ne crash pas 
        compute_properties("INVALIDE")

def test_compute_properties():
    result = compute_properties("CC(=O)Oc1ccccc1C(=O)O")
    assert result["log P"] == 1.3101
    assert result["Poids moléculaire"] == 180.15899999999996
    assert result["TPSA"] == 63.60000000000001
    assert result["Donneurs H"] == 1
    assert result["Accepteurs H"] == 3

def test_interpret_properties_insignifiant():
    result = interpret_properties({"log P": 1.0})
    assert result["log P"] == (1.0,"insignifiant")  

def test_interpret_properties_moderate():
    result = interpret_properties({"log P": 4.0})
    assert result["log P"] == (4.0,"moderate")

def test_interpret_properties_problematic():
    result = interpret_properties({"log P": 6.0})
    assert result["log P"] == (6.0,"problematic")

def test_properties_toxicity():
    input_dict = {
        "log P":             (1.0, "insignifiant"),
        "Poids moléculaire": (100, "insignifiant"),
        "TPSA":              (30,  "insignifiant"),
        "Donneurs H":        (1,   "insignifiant"),
        "Accepteurs H":      (2,   "insignifiant"),
    }
    result = properties_toxicity(input_dict)
    assert result == 0.0

