import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from utils import calculate_lipinski, extract_esol_descriptors

st.set_page_config(
    page_title="Drug Screening Dashboard",
    page_icon="🧪",
    layout="wide"
)

st.title("🧪 Drug Screening Dashboard")
st.markdown("""
This dashboard checks small molecules with **Lipinski's Rule of Five**
and pulls out the 5 features we use for solubility models.
""")

st.sidebar.header("Controls")
input_mode = st.sidebar.radio("Select Input Mode", ["Single SMILES", "Batch CSV Upload"])

if input_mode == "Single SMILES":
    st.sidebar.subheader("Single Molecule Input")
    default_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
    smiles_input = st.sidebar.text_input("Enter SMILES String:", value=default_smiles)

    if st.sidebar.button("Run Screening"):
        if not smiles_input:
            st.error("Please enter a SMILES string.")
        else:
            mol = Chem.MolFromSmiles(smiles_input)
            if mol is None:
                st.error("Invalid SMILES string. RDKit could not read it.")
            else:
                col1, col2 = st.columns([1, 2])

                with col1:
                    st.subheader("Molecule Picture")
                    img = Draw.MolToImage(mol, size=(300, 300))
                    st.image(img, caption="SMILES: " + smiles_input)

                with col2:
                    st.subheader("Lipinski Rule of Five")
                    lipinski_results = calculate_lipinski(smiles_input)
                    if lipinski_results:
                        metrics_df = pd.DataFrame(list(lipinski_results.items()), columns=["Parameter", "Value"])
                        metrics_df["Value"] = metrics_df["Value"].astype(str)
                        st.table(metrics_df)

                        if lipinski_results["Status"] == "Pass":
                            st.success("Status: Pass (drug-like)")
                        else:
                            st.warning("Status: Fail (not drug-like)")

                    st.subheader("Solubility Features")
                    descriptors = extract_esol_descriptors(smiles_input)
                    if descriptors:
                        desc_df = pd.DataFrame(list(descriptors.items()), columns=["Feature", "Value"])
                        desc_df["Value"] = desc_df["Value"].astype(str)
                        st.table(desc_df)

else:
    st.sidebar.subheader("Batch Upload")
    uploaded_file = st.sidebar.file_uploader("Upload a CSV file with a 'smiles' column", type=["csv"])

    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        if "smiles" not in df.columns:
            st.error("The CSV must have a column named 'smiles'.")
        else:
            st.info("Loaded " + str(len(df)) + " molecules.")
            if st.button("Run Batch Screening"):
                results = []
                for idx, row in df.iterrows():
                    smi = row["smiles"]
                    lip = calculate_lipinski(smi)
                    if lip:
                        lip["smiles"] = smi
                        results.append(lip)

                if results:
                    res_df = pd.DataFrame(results)
                    st.subheader("Batch Results")
                    st.dataframe(res_df)

                    csv_data = res_df.to_csv(index=False).encode("utf-8")
                    st.download_button(
                        label="Download Results (CSV)",
                        data=csv_data,
                        file_name="screening_report.csv",
                        mime="text/csv"
                    )
    else:
        st.info("Upload a CSV file with a 'smiles' column to start.")
