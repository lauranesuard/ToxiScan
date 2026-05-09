from rdkit import Chem
from rdkit.Chem import Descriptors
from src.toxiscan.toxicophores import find_toxicophores

def remove_redundant_toxicophores(detected_toxicophores: dict) -> dict:
    """
    Renvoie le dictionnaire de tous les groupes distincts trouvés dans une molécule.
    
    Parameters
    ----------
    found : dict
        les groupes toxicophores détectés par find_toxicophores
    
    Returns
    -------
    dict
        le nouveau dictionnaire sans les groupes detectés plusieurs fois par erreur
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
    Compte le nombre de toxicophores détectés dans une molécule.
    
    Parameters
    ----------
    smiles : str
        SMILES de la molécule à analyser.
    
    Returns
    -------
    int
        Nombre de toxicophores trouvés.
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
    Calcule la toxicité de la molécule
    
    Parameters
    ----------
    smiles : str
        SMILES de la molécule à analyser.
    
    Returns
    -------
    float
        Toxicité de la molécule par rapport à sa taille.
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
    Donne des informations complémeentaires sur la toxicité de la molécule
    
    Parameters
    ----------
    smiles : str
        SMILES de la molécule à analyser.
    
    Returns
    -------
    dict
        dictionnaire qui affichera 4 nouvelles propriétés pour chaque molécule
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