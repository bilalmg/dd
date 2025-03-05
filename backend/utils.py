from rdkit import Chem
from rdkit.Chem import Descriptors

def get_molecular_properties(smiles):
    """Calculate basic molecular properties from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {"error": "Invalid SMILES string"}

    properties = {
        "molecular_weight": Descriptors.MolWt(mol),
        "num_hydrogen_bond_donors": Descriptors.NumHDonors(mol),
        "num_hydrogen_bond_acceptors": Descriptors.NumHAcceptors(mol),
        "logP": Descriptors.MolLogP(mol),  # Lipophilicity
    }
    return properties
