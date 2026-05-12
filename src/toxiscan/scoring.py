from rdkit import Chem
from rdkit.Chem import Descriptors
from src.toxiscan.toxicophores import find_toxicophores

def remove_redundant_toxicophores(detected_toxicophores: dict) -> dict:
    """
    Returns a dictionary of all distinct groups found in a molecule.
    
    Parameters
    ----------
    found : dict
        the toxicophore groups detected by find_toxicophores
    
    Returns
    -------
    dict
        the new dictionary without the groups detected multiple times by mistake
    """
    redundant = []
    for name_A, atoms_A in detected_toxicophores.items():
        for name_B, atoms_B in detected_toxicophores.items():
            if name_A != name_B:
                if set(atoms_A) <= set(atoms_B):
                    redundant.append(name_A)

    for name in redundant:
        del detected_toxicophores[name]

    return detected_toxicophores

def count_toxicophores(smiles: str) -> int:
    """
    Counts the number of toxicophores detected in a molecule.
    
    Parameters
    ----------
    smiles : str
        SMILES string of the molecule to be analysed.
    
    Returns
    -------
    int
        Number of toxicophores found.
    """
    return len(remove_redundant_toxicophores(find_toxicophores(smiles)))

TOXICOPHORE_WEIGHTS = {
    # Nitrogen-based
    "Nitro group":       2,
    "Nitroso group":     2,
    "Primary aniline":   1,
    "Hydrazine":         2,
    "Aromatic azo":      2,
    "N-oxide":           1,
    "Hydroxamic acid":   1,
    # Carbonyl-based
    "Aldehyde":          2,
    "Acyl halide":       3,
    "Anhydride":         2,
    "Alpha halo ketone": 3,
    "Beta lactone":      3,
    # Electrophilic
    "Epoxide":           3,
    "Isocyanate":        3,
    "Michael acceptor":  2,
    "Aziridine":         3,
    "Activated alkyne":  2,
    # Halogen-based
    "Alkyl halide":      1,
    "Allylic halide":    2,
    # Sulfur-based
    "Thiol":             1,
    "Thiocarbonyl":      1,
    # Peroxide
    "Peroxide":          3,
    # Polycyclic aromatic
    "Quinone":           2,
    "Coumarin":          1,
    # Reactive oxygen
    "Hydroperoxide":     3,
}

# les poids sont basés sur la réactivité électrophile connue des groupes fonctionnels

def toxicity_approximation(smiles: str) -> float:
    """
    Calculates the toxicity of the molecule
    
    Parameters
    ----------
    smiles : str
        SMILES string for the molecule to be analysed.
    
    Returns
    -------
    float
        Toxicity of the molecule relative to its size.
    """
    clean_toxicophores = remove_redundant_toxicophores(find_toxicophores(smiles))
    score = 0

    for name, atoms in clean_toxicophores.items(): 
        score += TOXICOPHORE_WEIGHTS[name]
    
    number_atoms = Chem.MolFromSmiles(smiles).GetNumAtoms()
    score = score / number_atoms 
    
    return score

def compute_properties(smiles: str) -> dict:
    """
    Provides additional information on the toxicity of the molecule
    
    Parameters
    ----------
    smiles : str
        SMILES string for the molecule to be analysed.
    
    Returns
    -------
    dict
        dictionary containing 4 new properties for each molecule
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES invalide : '{smiles}'")
    
    return {
        "log P": Descriptors.MolLogP(mol),
        "Poids moléculaire": Descriptors.MolWt(mol),
        "TPSA": Descriptors.TPSA(mol),
        "Donneurs H": Descriptors.NumHDonors(mol),
        "Accepteurs H": Descriptors.NumHAcceptors(mol),
    }

THRESHOLDS = {
    "log P":              (3, 5),
    "Poids moléculaire":  (300, 500),
    "TPSA":               (60, 140),
    "Donneurs H":         (2, 5),
    "Accepteurs H":       (5, 10),
}

def interpret_properties(properties: dict) -> dict:
    """
    Interprets the physicochemical properties of a molecule
    and indicates whether each value falls within a favourable, 
    moderate or problematic range according to Lipinski’s rules.

    Parameters
    ----------
    properties : dict
        dictionary returned by compute_properties(), containing
        the values for logP, molecular weight, TPSA, H donors
        and H acceptors.

    Returns
    -------
    dict
        dictionary with an indication for each property:
        “insignificant”, “moderate” or “problematic”.
    """

    scale = {}

    for name, values in properties.items():
        low, high = THRESHOLDS[name]
        if values <= low : 
            scale[name] = (values, "insignifiant")
        
        elif low < values <= high:
            scale[name] = (values,"moderate")
        
        else:
            scale[name] = (values,"problematic")
        
    return scale

