from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import display, Image
import py3Dmol

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

def draw_molecule_3d(smiles: str, toxicophores_found: dict) -> str:
    """
    Generate an interactive 3D visualization of a molecule with toxic
    fragments highlighted in green.
    
    Parameters :
        smiles : str
            the SMILES string of the molecule
        toxicophores_found : dict
            dictionary returned by find_toxicophores(), with toxicophore
            names as keys and atom indices as values
    
    Returns :
        str
            HTML string of the 3D viewer
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: '{smiles}'")
    
    # Generate 3D coordinates
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)

    mol = Chem.RemoveHs(mol) 
    mol_block = Chem.MolToMolBlock(mol)
    
    # Collect toxic atom indices
    highlight_atoms = []
    for indices in toxicophores_found.values():
        highlight_atoms.extend(indices)
    highlight_atoms = list(set(highlight_atoms))
    
    # Create viewer
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(mol_block, "mol")
    
    # Style all atoms as ballstick
    viewer.setStyle({"stick": {"radius": 0.15}, "sphere": {"scale": 0.3}})
    
    # Highlight toxic atoms in green
    for idx in highlight_atoms:
        viewer.setStyle({"serial": idx}, 
                       {"stick": {"radius": 0.15, "color": "#7FFF00"},
                        "sphere": {"scale": 0.4, "color": "#7FFF00"}})
    
    viewer.setBackgroundColor("white")
    viewer.zoomTo()
    
    return viewer._make_html()
