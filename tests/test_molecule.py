from unittest.mock import patch, MagicMock
from src.toxiscan.molecule import get_smiles, get_molecule_info

def test_gets_smiles_found(): 
    mock_response = MagicMock()    # création du faux objet, MagicMock s'applique à ce qu'on met dans notre fonction
    mock_response.status_code = 200 # fais semblant d'avoir réussi
    mock_response.json.return_value = {
        "PropertyTable": {
            "Properties": [{"ConnectivitySMILES": "CC(=O)Oc1ccccc1C(=O)O"}]
        }
    }
    with patch("src.toxiscan.molecule.requests.get", return_value = mock_response):
        result = get_smiles("aspirin")
        assert result == "CC(=O)Oc1ccccc1C(=O)O"

def test_get_smiles_not_found():
    mock_response = MagicMock()
    mock_response.status_code = 404
    with patch("src.toxiscan.molecule.requests.get", return_value=mock_response):
        result = get_smiles("moleculeinventee")
        assert "not found" in result

def test_get_molecule_info_found():
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {
        "PropertyTable": {
            "Properties": [{
                "ConnectivitySMILES": "CC(=O)Oc1ccccc1C(=O)O",
                "MolecularFormula": "C9H8O4",
                "MolecularWeight": "180.16",
                "IUPACName": "2-acetyloxybenzoic acid",
            }]
        }
    }
    with patch("src.toxiscan.molecule.requests.get", return_value=mock_response):
        result = get_molecule_info("aspirin")
        assert result["name"] == "aspirin"
        assert result["formula"] == "C9H8O4"
        assert result["molecular_weight"] == "180.16"
        assert result["iupac_name"] == "2-acetyloxybenzoic acid"
        assert result["smiles"] == "CC(=O)Oc1ccccc1C(=O)O"

def test_get_molecule_info_not_found():
    mock_response = MagicMock()
    mock_response.status_code = 404
    with patch("src.toxiscan.molecule.requests.get", return_value=mock_response):
        result = get_molecule_info("moleculeinventee")
        assert "error" in result 
