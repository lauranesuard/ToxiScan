import streamlit as st
import sys
sys.path.insert(0, '.')
from src.toxiscan.molecule import get_smiles
from src.toxiscan.toxicophores import find_toxicophores
from src.toxiscan.scoring import remove_redundant_toxicophores
from src.toxiscan.visualization import draw_molecule, draw_molecule_3d

st.set_page_config(page_title="ToxiScan", layout="wide")
st.title("🧪 ToxiScan — Molecular Toxicity Analyzer")

molecule_name = st.text_input("Enter a molecule name:", placeholder="e.g. aspirin, caffeine, ethanol")

if molecule_name:
    smiles = get_smiles(molecule_name)
    
    if "not found" in str(smiles):
        st.error(smiles)
    else:
        st.success(f"SMILES: {smiles}")
        
        col1, col2 = st.columns(2)
        
        # 2D visualization
        with col1:
            st.subheader("2D Structure")
            tox = find_toxicophores(smiles)
            tox_clean = remove_redundant_toxicophores(tox)
            draw_molecule(smiles, tox_clean)
            st.image("molecule.svg")
        
        # 3D visualization
        with col2:
            st.subheader("3D Structure")
            html_3d = draw_molecule_3d(smiles, tox_clean)
            st.components.v1.html(html_3d, height=400)