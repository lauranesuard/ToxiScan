from rdkit import Chem
TOXICOPHORES = {
    # Nitrogen-based
    "Nitro group":      "[N+](=O)[O-]",
    "Nitroso group":    "[N]=O",
    "Primary aniline":  "[NH2]c1ccccc1",
    "Hydrazine":        "[NH]-[NH2]",

    # Carbonyl-based
    "Aldehyde":         "[CH]=O",
    "Acyl halide":      "C(=O)[F,Cl,Br,I]",
    "Anhydride":        "C(=O)OC(=O)",

    # Electrophilic
    "Epoxide":          "C1OC1",
    "Isocyanate":       "N=C=O",
    "Michael acceptor": "C=CC=O",

    # Halogen-based
    "Alkyl halide":     "[CX4][F,Cl,Br,I]",
    "Allylic halide":   "C=C[CH2][F,Cl,Br,I]",

    # Sulfur-based
    "Thiol":            "[SH]",
    "Thiocarbonyl":     "C=S",

    # Peroxide
    "Peroxide":         "OO",
}

def find_toxicophores(smiles: str) -> dict :
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

    

