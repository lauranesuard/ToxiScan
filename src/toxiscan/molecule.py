import requests

def get_smiles(molecule_name: str) -> str:
    """
    Convert a molecule name to its SMILES representation using PubChem API.
    
    Parameters : 
        molecule_name : str
            the name of the molecule (e.g. "aspirin", "caffeine")
    
    Returns : 
        str : the SMILES string of the molecule, or None if not found
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{molecule_name}/property/CanonicalSMILES/JSON"
    
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        smiles = data["PropertyTable"]["Properties"][0]["ConnectivitySMILES"]
        return smiles
    else:
        return f"Molecule '{molecule_name}' not found in PubChem. Please check the name and try again."
    
