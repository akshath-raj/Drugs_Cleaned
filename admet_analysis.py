import os
import glob
import pandas as pd
from rdkit import Chem
from adme_py import ADME
from admet_ai import ADMETModel

# Define the columns to display in the final UI
DISPLAY_COLUMNS = [
    "Filename",
    # Physicochemical (from adme-py)
    "MW", "WLogP", "TPSA", "HBD", "HBA", "Rotatable Bonds", "Fsp3",
    # Absorption / Distribution
    "GI Absorption", "Caco-2 (Wang)", "BBB (Martins)", "PPB (AZ)",
    # Metabolism
    "CYP3A4 Inhibition", "CYP2D6 Inhibition",
    # Toxicity
    "hERG", "Ames", "DILI", "Carcinogenicity", "LD50 (Zhu)",
    # Developability
    "Lipinski", "PAINS", "Brenk", "QED", "SA Score"
]

def pdb_to_smiles(pdb_path):
    """Convert PDB file to SMILES string using RDKit."""
    try:
        mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol)
    except Exception:
        return None

def run_admet_prediction():
    """
    Scans docking_results/pdb for *ligand_poses.pdb files,
    runs ADME (SwissADME style) and ADMET-AI analysis, 
    and returns a filtered DataFrame and CSV path.
    """
    
    # 1. Find Files
    search_path = os.path.join("docking_results", "pdb", "*ligand_poses.pdb")
    pdb_files = glob.glob(search_path)
    
    if not pdb_files:
        return "No ligand PDB files found in docking_results/pdb. Please run docking first.", None, None

    # 2. Convert to SMILES
    compounds = []
    smiles_list = []
    
    for pdb_file in pdb_files:
        smiles = pdb_to_smiles(pdb_file)
        if smiles:
            filename = os.path.basename(pdb_file)
            compounds.append({"Filename": filename, "SMILES": smiles})
            smiles_list.append(smiles)

    if not compounds:
        return "Could not generate valid SMILES from PDB files.", None, None

    # 3. Run adme-py (Physicochemical, SA Score, etc.)
    adme_data = []
    
    for cmp in compounds:
        row = {}
        try:
            res = ADME(cmp["SMILES"]).calculate()
            
            # Physiochemical
            physio = res.get("physiochemical", {})
            row["MW"] = physio.get("molecular_weight", "N/A")
            row["TPSA"] = physio.get("tpsa", "N/A")
            row["HBD"] = physio.get("num_h_donors", "N/A")
            row["HBA"] = physio.get("num_h_acceptors", "N/A")
            row["Rotatable Bonds"] = physio.get("num_rotatable_bonds", "N/A")
            row["Fsp3"] = physio.get("sp3_carbon_ratio", "N/A")
            
            # Lipophilicity
            row["WLogP"] = res.get("lipophilicity", {}).get("wlogp", "N/A")
            
            # Pharmacokinetics
            row["GI Absorption"] = res.get("pharmacokinetics", {}).get("gastrointestinal_absorption", "N/A")
            
            # Druglikeness / Medicinal
            row["Lipinski"] = res.get("druglikeness", {}).get("lipinski", "N/A")
            
            # Handle Booleans for PAINS/Brenk
            pains = res.get("medicinal", {}).get("pains", False)
            row["PAINS"] = "Yes" if pains else "No"
            
            brenk = res.get("medicinal", {}).get("brenk", False)
            row["Brenk"] = "Yes" if brenk else "No"
            
            # SA Score (The one you asked about)
            row["SA Score"] = res.get("medicinal", {}).get("synthetic_accessibility", "N/A")
            
        except Exception as e:
            print(f"ADME error for {cmp['Filename']}: {e}")
            row.update({k: "Error" for k in ["MW", "WLogP", "TPSA", "Lipinski", "SA Score"]})

        adme_data.append(row)

    df_adme = pd.DataFrame(adme_data)

    # 4. Run ADMET-AI (Toxicity, Metabolism, QED)
    try:
        model = ADMETModel()
        preds_df = model.predict(smiles_list)
        
        admet_mapped = pd.DataFrame()
        
        # Mapping columns exactly from your requested list
        admet_mapped["Caco-2 (Wang)"] = preds_df.get("Caco2_Wang", "N/A")
        admet_mapped["BBB (Martins)"] = preds_df.get("BBB_Martins", "N/A")
        admet_mapped["PPB (AZ)"] = preds_df.get("PPBR_AZ", "N/A")
        
        admet_mapped["CYP3A4 Inhibition"] = preds_df.get("CYP3A4_Veith", "N/A")
        admet_mapped["CYP2D6 Inhibition"] = preds_df.get("CYP2D6_Veith", "N/A")
        
        admet_mapped["hERG"] = preds_df.get("hERG", "N/A")
        admet_mapped["Ames"] = preds_df.get("AMES", "N/A")
        admet_mapped["DILI"] = preds_df.get("DILI", "N/A")
        admet_mapped["Carcinogenicity"] = preds_df.get("Carcinogens_Lagunin", "N/A")
        admet_mapped["LD50 (Zhu)"] = preds_df.get("LD50_Zhu", "N/A")
        
        admet_mapped["QED"] = preds_df.get("QED", "N/A")

    except Exception as e:
        print(f"ADMET-AI Error: {e}")
        admet_mapped = pd.DataFrame(index=range(len(compounds)))

    # 5. Merge Data
    df_base = pd.DataFrame(compounds)
    
    # Concatenate (reset index to ensure alignment)
    final_df = pd.concat([
        df_base.reset_index(drop=True),
        df_adme.reset_index(drop=True),
        admet_mapped.reset_index(drop=True)
    ], axis=1)

    # 6. Filter and Reorder Columns
    existing_cols = [c for c in DISPLAY_COLUMNS if c in final_df.columns]
    final_df_clean = final_df[existing_cols]

    # Round numeric columns for cleaner UI
    numeric_cols = ["MW", "WLogP", "TPSA", "Fsp3", "SA Score", "QED", "LD50 (Zhu)", "PPB (AZ)"]
    for col in numeric_cols:
        if col in final_df_clean.columns:
            final_df_clean[col] = pd.to_numeric(final_df_clean[col], errors='coerce').round(2)

    # 7. Save to CSV
    os.makedirs("results", exist_ok=True)
    csv_path = os.path.join("results", "final_admet_report.csv")
    final_df_clean.to_csv(csv_path, index=False)

    return f"Analysis Complete. Processed {len(final_df_clean)} ligands.", final_df_clean, csv_path