from toxiscan.scoring import remove_redundant_toxicophores, count_toxicophores, compute_properties, interpret_properties, properties_toxicity, toxicity_approximation
import pytest

def test_remove_redundant_no_redundancy():
    """
    Test that remove_redundant_toxicophores keeps all groups
    when there is no redundancy between them.
    """
    input_dict = {
        "Epoxide": (1, 2, 3),
        "Nitro group": (4, 5, 6),
    }
    result = remove_redundant_toxicophores(input_dict)
    assert "Epoxide" in result
    assert "Nitro group" in result

def test_remove_redundant_with_redundancy():
    """
    Test that remove_redundant_toxicophores removes a group
    whose atoms are a subset of another group's atoms.
    """
    input_dict = {
        "Epoxide": (1, 2),
        "Nitro group": (1, 2, 3, 4),
    }
    result = remove_redundant_toxicophores(input_dict)
    assert "Nitro group" in result
    assert "Epoxide" not in result

def test_count_toxicophores_zero():
    """
    Test that count_toxicophores returns 0 for aspirin,
    which contains no known toxicophores.
    """
    result = count_toxicophores("CC(=O)Oc1ccccc1C(=O)O")  # aspirine
    assert result == 0

def test_count_toxicophores_nonzero():
    """
    Test that count_toxicophores returns 1 for formaldehyde,
    which contains one aldehyde group.
    """
    result = count_toxicophores("C=O")  # formaldéhyde, 1 aldéhyde
    assert result == 1

def test_toxicity_approximation_zero():
    """
    Test that toxicity_approximation returns 0.0 for aspirin,
    which contains no toxicophores.
    """
    result = toxicity_approximation("CC(=O)Oc1ccccc1C(=O)O")  
    assert result == 0.0

# on re test la fonction dans le cas où un toxicophore est détecté
def test_toxicity_approximation_nonzero():
    """
    Test that toxicity_approximation returns a non-zero score
    for formaldehyde, which contains an aldehyde group.
    """
    result = toxicity_approximation("C=O")  # formaldéhyde, aldéhyde
    assert result == 1.0

def test_compute_properties_fail():
    """
    Test that compute_properties raises a ValueError
    when given an invalid SMILES string.
    """
    with pytest.raises(ValueError): # pour signifier je m'attends à une ValueError -> normal, le code ne crash pas 
        compute_properties("INVALIDE")

def test_compute_properties():
    """
    Test that compute_properties returns correct physicochemical
    properties for aspirin.
    """
    result = compute_properties("CC(=O)Oc1ccccc1C(=O)O")
    assert result["log P"] == 1.3101
    assert result["Molecular weigh"] == 180.15899999999996
    assert result["TPSA"] == 63.60000000000001
    assert result["H Donors"] == 1
    assert result["H Acceptors"] == 3

def test_interpret_properties_insignifiant():
    """
    Test that interpret_properties returns 'insignifiant'
    for a logP value below the lower threshold.
    """
    result = interpret_properties({"log P": 1.0})
    assert result["log P"] == (1.0,"insignifiant")  

def test_interpret_properties_moderate():
    """
    Test that interpret_properties returns 'moderate'
    for a logP value between the two thresholds.
    """
    result = interpret_properties({"log P": 4.0})
    assert result["log P"] == (4.0,"moderate")

def test_interpret_properties_problematic():
    """
    Test that interpret_properties returns 'problematic'
    for a logP value above the upper threshold.
    """
    result = interpret_properties({"log P": 6.0})
    assert result["log P"] == (6.0,"problematic")

def test_properties_toxicity():
    """
    Test that properties_toxicity returns 0.0 when all
    properties are classified as 'insignifiant'.
    """
    input_dict = {
        "log P":             (1.0, "insignifiant"),
        "Molecular weigh": (100, "insignifiant"),
        "TPSA":              (30,  "insignifiant"),
        "H Donors":        (1,   "insignifiant"),
        "H Acceptors":      (2,   "insignifiant"),
    }
    result = properties_toxicity(input_dict)
    assert result == 0.0

