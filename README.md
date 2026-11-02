# Drug Screening Dashboard

A small web app that checks molecules with Lipinski's Rule of Five and pulls out features used for solubility models.

## What it does
- Single mode: type a SMILES, see a picture of the molecule, its Lipinski results, and 5 solubility features (weight, H-bond donors, rotatable bonds, rings, TPSA).
- Batch mode: upload a CSV with a smiles column, screen many molecules at once, and download the results.
- Bad SMILES strings are skipped instead of crashing the app.

## Tools used
Python, RDKit, Streamlit, pandas

## How to run it
pip3 install streamlit pandas rdkit pillow
python3 -m streamlit run app.py

## Files
- app.py: the web page
- utils.py: the Lipinski and feature functions
- test_molecules.csv: a small test file (one fake molecule on purpose)

## Notes
- Pass means 0 or 1 Lipinski violations, as in the original rule.
- Built while learning computational chemistry. Part of my portfolio with a solubility predictor: github.com/nkosanabrendon-bit/solubility-predictor
