# ToxiScan 🧪
[logo/banner]
[badges]

## Package description 📖
**ToxiScan** is a Python package designed to evaluate the toxicity of molecules, aimed at chemists, students, and researchers who need fast and visual toxicity insights.

Given a molecule name as input, ToxiScan:
  1. Retrieves its SMILES representation from the **PubChem** database
  2. Detects toxic functional groups via **SMARTS pattern matching**
  3. Visualizes the molecule in **2D and 3D** with toxic fragments highlighted in green
  4. Computes a **toxicity score** based on detected toxicophores and physicochemical properties

## Authors 🧑‍🎓
This project was developed as part of the *Practical Programming in Chemistry* course at EPFL.

| Name | GitHub |
|------|--------|
| Laurane Suard | [@lauranesuard](https://github.com/lauranesuard) |
| Flore Leveillé-Nizerolle | [@floreln](https://github.com/floreln) |
| Nina Deruaz | [@ninaderuaz](https://github.com/ninaderuaz) |

## Features 🔬

- **Toxicophore detection**: scan of the molecule (SMILES string form) to identify the toxic substructures (epoxides, nitro groups, aldehydes...)
- **Toxicity scoring**: normalized score based on type and severity of detected toxicophores, combined with physicochemical properties (Lipinski's rules)
- **2D visualization**: molecule drawn in black and white with toxic atoms highlighted in green
- **3D interactive visualization**: rotate and zoom on the molecule with py3Dmol, toxic atoms highlighted in green
- **Streamlit app**: user-friendly interface to analyze any molecule in seconds

## Installation 💻
First, clone the repository and navigate into the project directory:

```bash
git clone https://github.com/lauranesuard/ToxiScan.git
cd ToxiScan
```

Then, create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate toxiscan-env
```

Finally, install the package in editable mode:

```bash
pip install -e .
```

## Requirements 📝

## Usage 🚀

## Interface 🌐
[screenshots Streamlit]

## Run tests ✅

## Troubleshooting 🔧

## License 📜

