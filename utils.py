from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors


def calculate_lipinski(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    weight = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    donors = Lipinski.NumHDonors(mol)
    acceptors = Lipinski.NumHAcceptors(mol)

    violations = 0
    if weight >= 500:
        violations += 1
    if logp >= 5:
        violations += 1
    if donors >= 5:
        violations += 1
    if acceptors >= 10:
        violations += 1

    if violations <= 1:
        status = "Pass"
    else:
        status = "Fail"

    return {
        "Molecular Weight": round(weight, 2),
        "LogP": round(logp, 2),
        "H-Bond Donors": donors,
        "H-Bond Acceptors": acceptors,
        "Violations": violations,
        "Status": status,
    }


def extract_esol_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return {
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "H-Bond Donors": Lipinski.NumHDonors(mol),
        "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),
        "Rings": rdMolDescriptors.CalcNumRings(mol),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 2),
    }
