import streamlit as st
import sys
sys.path.insert(0, '.')
import pandas as pd
from src.toxiscan.scoring import toxicity_approximation, compute_properties, interpret_properties, TOXICOPHORE_WEIGHTS
from src.toxiscan.molecule import get_smiles
from src.toxiscan.toxicophores import find_toxicophores
from src.toxiscan.scoring import remove_redundant_toxicophores
from src.toxiscan.visualization import draw_molecule, draw_molecule_3d

st.set_page_config(page_title="ToxiScan", layout="wide")
st.title("🧪 ToxiScan — Molecular Toxicity Analyzer")

molecule_name = st.text_input("Enter a molecule name:", placeholder="e.g. aspirin, caffeine, ethanol", key="molecule_input")

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
        
        st.divider()
        st.subheader("🔬 Detected Toxic Groups")

        if not tox_clean:
            st.success("✅ No toxic substructures detected.")
        else:
            rows = []
            for name in tox_clean:
                rows.append({"Group": name, "Weight": TOXICOPHORE_WEIGHTS[name]})
            
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

        st.divider()
        st.subheader("☠️ Toxicity Score")

        score = toxicity_approximation(smiles)

        if score == 0:
            level, icon = "Non-toxic", "🟢"
        elif score < 0.1:
            level, icon = "Low toxicity", "🟡"
        elif score < 0.25:
            level, icon = "Moderate toxicity", "🟠"
        else:
            level, icon = "High toxicity", "🔴"

        st.metric("Score", f"{score:.3f}")
        st.progress(min(score / 0.5, 1.0), text=f"{icon} {level}")

        st.divider()
        st.subheader("📊 Physicochemical Properties")
        st.caption ("These properties complement the toxicity score: while the score is based on detected toxic groups, these physicochemical parameters (based on Lipinski's rules) give a broader picture of the molecule's behavior.")

        props = compute_properties(smiles)
        interpreted = interpret_properties(props)

        cols = st.columns(5)
        for col, (name, (value, level)) in zip(cols, interpreted.items()):
            col.metric(label=name, value=round(value, 2), delta=level, delta_color="off")
