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
    

def get_molecule_info(molecule_name: str) -> dict:
    """
    Retrieve basic information about a molecule from PubChem.
    
    Parameters : 
        molecule_name : str
            the name of the molecule (e.g. "aspirin", "caffeine")
    
    Returns :
        dict
            A dictionary containing basic molecular information,
            or an error message if the molecule is not found.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{molecule_name}/property/ConnectivitySMILES,MolecularFormula,MolecularWeight,IUPACName/JSON"

    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        props = data["PropertyTable"]["Properties"][0]
        return {
            "name": molecule_name,
            "formula": props["MolecularFormula"],
            "molecular_weight": props["MolecularWeight"],
            "iupac_name": props["IUPACName"],
            "smiles": props["ConnectivitySMILES"]
        }
    else:
        return {"error": f"Molecule '{molecule_name}' not found in PubChem. Please check the name and try again."}
    
