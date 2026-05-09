from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import display, Image

def draw_molecule(smiles: str, toxicophores_found: dict) -> None:
    """
    Draw a molecule with toxic fragments highlighted in yellow-green.
    
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

    # Draw with black atoms and yellow-green highlights
    drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
    drawer.drawOptions().addAtomIndices = False
    drawer.drawOptions().useBWAtomPalette()
    
    highlight_color = {atom: (0.6, 0.9, 0.2) for atom in highlight_atoms}
    
    drawer.DrawMolecule(mol, 
                        highlightAtoms=highlight_atoms,
                        highlightAtomColors=highlight_color,
                        highlightBonds=[])
    drawer.FinishDrawing()
    
    svg = drawer.GetDrawingText()
    with open("molecule.svg", "w") as f:
        f.write(svg)