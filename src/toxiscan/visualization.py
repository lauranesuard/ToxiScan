from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import display, Image

def draw_molecule(smiles: str, toxicophores_found: dict) -> None:
    """
    Draw a molecule with toxic fragments highlighted in red.
    
    Parameters : 
        smiles : str
            the SMILES string of the molecule
        toxicophores_found : dict
            dictionary returned by find_toxicophores(), with toxicophore
            names as keys and atom indices as values
    
    Returns : 
        None
            Displays the molecule image inline
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: '{smiles}'")
    
    # Collect all toxic atom indices
    highlight_atoms = []
    for indices in toxicophores_found.values():
        highlight_atoms.extend(indices)
    highlight_atoms = list(set(highlight_atoms))
    
    # Draw molecule
    img = Draw.MolToImage(mol, size=(400, 300), highlightAtoms=highlight_atoms)
    img.save("molecule.png")
    display(Image("molecule.png"))