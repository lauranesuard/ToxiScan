from rdkit import Chem

TOXICOPHORES = {
    # Nitrogen-based
    "Nitro group":      "[N+](=O)[O-]",
    "Nitroso group":    "[N]=O",
    "Primary aniline":  "[NH2]c1ccccc1",
    "Hydrazine":        "[NH]-[NH2]",
    "Aromatic azo":     "c-N=N-c",
    "N-oxide":          "[N+][O-]",
    "Hydroxamic acid":  "C(=O)NO",

    # Carbonyl-based
    "Aldehyde":          "[CX3;H1,H2](=O)",
    "Acyl halide":       "C(=O)[F,Cl,Br,I]",
    "Anhydride":         "C(=O)OC(=O)",
    "Alpha halo ketone": "C(=O)C[F,Cl,Br,I]",
    "Beta lactone":      "C1CC(=O)O1",

    # Electrophilic
    "Epoxide":          "C1OC1",
    "Isocyanate":       "N=C=O",
    "Michael acceptor": "C=CC=O",
    "Aziridine":        "C1NC1",
    "Activated alkyne": "C#CC=O",

    # Halogen-based
    "Alkyl halide":     "[CX4][F,Cl,Br,I]",
    "Allylic halide":   "C=C[CH2][F,Cl,Br,I]",

    # Sulfur-based
    "Thiol":            "[SH]",
    "Thiocarbonyl":     "C=S",

    # Peroxide
    "Peroxide":         "OO",

    #Polycyclic aromatic
    "Quinone": "O=C1C=CC(=O)C=C1",
    "Coumarin": "O=C1OC2=CC=CC=C2C=C1",

    #Reactive oxygen
    "Hydroperoxide": "[OX2][OX2H]",
}

def find_toxicophores(smiles: str) -> dict:
    """
    Detect toxic functional groups (toxicophores) in a molecule.
    
    Parameters : 
        smiles : str
            the SMILES string of the molecule
    
    Returns : 
        dict
            a dictionary where keys are toxicophore names and values are 
            lists of atom indices where the toxicophore was found.
            Returns an empty dict if no toxicophores are detected.
    
    Raises
    ------
    ValueError
        If the SMILES string is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES invalide : '{smiles}'")
    
    results = {}

    for name, smarts in TOXICOPHORES.items():
        pattern = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pattern)

        if matches:
            atom_indices = list(set(idx for match in matches for idx in match))
            results[name] = atom_indices
            
    return results

    

