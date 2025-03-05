import streamlit as st
import requests

# Icons
icon_molecule = "🔬"
icon_bioactivity = "🧪"
icon_interactions = "⚗️"
icon_warning = "⚠️"

# Streamlit UI
st.title("🧬 Drug Discovery Assistant")

# User input
st.image("https://upload.wikimedia.org/wikipedia/commons/3/3f/Chemical_DNA.png", width=100)
smiles = st.text_input("Enter SMILES notation for the molecule:")

if st.button("Analyze"):
    if smiles:
        response = requests.post("http://127.0.0.1:5000/process_molecule", json={"molecule": smiles})

        if response.status_code == 200:
            result = response.json()

            st.subheader("🔍 Analysis Results:")

            # Check for errors before displaying results
            if "Molecular Properties" in result and result["Molecular Properties"] != "Invalid molecule":
                st.subheader(f"{icon_molecule} Molecular Properties")
                st.json(result["Molecular Properties"])
            else:
                st.error(f"{icon_warning} Invalid molecule. Please check your SMILES input.")

            if "AI Bioactivity Prediction" in result:
                st.subheader(f"{icon_bioactivity} AI Bioactivity Prediction")
                st.write(result["AI Bioactivity Prediction"])

            if "Known Drug Interactions" in result:
                st.subheader(f"{icon_interactions} Known Drug Interactions")
                st.write(", ".join(result["Known Drug Interactions"]))
        else:
            st.error(f"{icon_warning} Error: Could not process the request. Check the backend.")
    else:
        st.error(f"{icon_warning} Please enter a valid SMILES notation.")
