from flask import Flask, request, jsonify
import os
import openai
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors

app = Flask(__name__)

# Load OpenAI API Key from environment variable
openai.api_key = "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"

# ChEMBL clients
drug_interaction = new_client.molecule

# Function to calculate molecular properties
def get_molecular_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    properties = {
        "Molecular Weight": round(Descriptors.MolWt(mol), 3),
        "H-Bond Acceptors": Descriptors.NumHAcceptors(mol),
        "H-Bond Donors": Descriptors.NumHDonors(mol),
        "logP": round(Descriptors.MolLogP(mol), 3)
    }
    return properties

# Function to predict bioactivity using OpenAI
def predict_bioactivity(smiles):
    try:
        response = openai.chat.completions.create(
            model="gpt-4",
            messages=[{"role": "user", "content": f"Predict the bioactivity of the molecule with SMILES: {smiles}"}]
        )
        return response.choices[0].message.content.strip()
    except Exception as e:
        return f"Error in AI prediction: {str(e)}"

# Function to check known drug interactions
def get_drug_interactions(smiles):
    try:
        results = drug_interaction.filter(molecule_structures__canonical_smiles=smiles)
        interactions = [res["pref_name"] for res in results if "pref_name" in res]
        return interactions if interactions else ["No known interactions found"]
    except Exception as e:
        return [f"Error fetching interactions: {str(e)}"]

@app.route('/process_molecule', methods=['POST'])
def process_molecule():
    data = request.json
    smiles = data.get("molecule")
    
    if not smiles:
        return jsonify({"error": "No molecule provided"}), 400
    
    properties = get_molecular_properties(smiles)
    if properties is None:
        return jsonify({"error": "Invalid molecule"}), 400
    
    bioactivity = predict_bioactivity(smiles)
    interactions = get_drug_interactions(smiles)
    
    response = {
        "Molecular Properties": properties,
        "AI Bioactivity Prediction": bioactivity,
        "Known Drug Interactions": interactions
    }
    return jsonify(response)

if __name__ == '__main__':
    app.run(debug=True)
