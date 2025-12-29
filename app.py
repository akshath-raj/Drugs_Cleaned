"""
Main Gradio interface for Protein Structure Finder & Analyzer
"""

import os
import requests
import tempfile
import gradio as gr
import pandas as pd

# Import modules
from config import current_pdb_info, PROTEINS_DIR
from utils import map_disease_to_protein, search_pdb_for_first_hit, remove_ligands_from_pdb
from visualization import show_structure
from ramachandran import run_ramplot
from prankweb import run_prankweb_prediction
from protein_prep import prepare_protein_meeko
from docking import run_molecular_docking, display_docked_structure
# Import the updated ADMET module
from admet_analysis import run_admet_prediction


def process_disease(disease_name: str):
    """Main function to process disease and return structure."""
    
    if not disease_name.strip():
        current_pdb_info.update({"pdb_id": None, "pdb_path": None, "prepared_pdbqt": None, "docking_results": None, "prankweb_csv": None})
        return {
            info_box: gr.update(visible=False),
            structure_viewer: gr.update(value=""),
            download_file: gr.update(value=None),
            search_status: gr.update(value="⚠️ Please enter a disease or condition", visible=True)
        }
    
    # Map disease to protein
    protein_name = map_disease_to_protein(disease_name)
    
    if not protein_name:
        current_pdb_info.update({"pdb_id": None, "pdb_path": None, "prepared_pdbqt": None, "docking_results": None, "prankweb_csv": None})
        return {
            info_box: gr.update(visible=False),
            structure_viewer: gr.update(value=""),
            download_file: gr.update(value=None),
            search_status: gr.update(value="❌ No protein mapping found", visible=True)
        }
    
    # Search PDB
    pdb_id = search_pdb_for_first_hit(protein_name)
    
    if not pdb_id:
        current_pdb_info.update({"pdb_id": None, "pdb_path": None, "prepared_pdbqt": None, "docking_results": None, "prankweb_csv": None})
        return {
            info_box: gr.update(visible=False),
            structure_viewer: gr.update(value=""),
            download_file: gr.update(value=None),
            search_status: gr.update(value="❌ No PDB structure found", visible=True)
        }
    
    # Download PDB file
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        pdb_content = response.text
        
        # Clean structure
        pdb_content, stats = remove_ligands_from_pdb(pdb_content, 'A')
        
        # Save to proteins folder
        proteins_folder = PROTEINS_DIR
        os.makedirs(proteins_folder, exist_ok=True)
        pdb_path = os.path.join(proteins_folder, f"{pdb_id}.pdb")
        
        with open(pdb_path, 'w') as f:
            f.write(pdb_content)
        
        # Update global variable
        current_pdb_info.update({"pdb_id": pdb_id, "pdb_path": pdb_path, "prepared_pdbqt": None, "docking_results": None, "prankweb_csv": None})
        
        # Build info display
        info_html = f"""
        <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 24px; border-radius: 16px; color: white; box-shadow: 0 8px 32px rgba(0,0,0,0.1);">
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px;">
                <div>
                    <div style="font-size: 13px; opacity: 0.9; font-weight: 600; text-transform: uppercase; letter-spacing: 0.5px; margin-bottom: 8px;">Disease/Condition</div>
                    <div style="font-size: 20px; font-weight: 700;">{disease_name}</div>
                </div>
                <div>
                    <div style="font-size: 13px; opacity: 0.9; font-weight: 600; text-transform: uppercase; letter-spacing: 0.5px; margin-bottom: 8px;">Target Protein</div>
                    <div style="font-size: 20px; font-weight: 700;">{protein_name}</div>
                </div>
                <div>
                    <div style="font-size: 13px; opacity: 0.9; font-weight: 600; text-transform: uppercase; letter-spacing: 0.5px; margin-bottom: 8px;">PDB Structure ID</div>
                    <div style="font-size: 20px; font-weight: 700;">{pdb_id}</div>
                </div>
            </div>
        </div>
        """
        
        # Create 3D visualization
        structure_html = show_structure(pdb_content, pdb_id, protein_name)
        
        # Create download file
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
        temp_file.write(pdb_content)
        temp_file.close()
        
        return {
            info_box: gr.update(value=info_html, visible=True),
            structure_viewer: gr.update(value=structure_html),
            download_file: gr.update(value=temp_file.name),
            search_status: gr.update(value="✅ Structure loaded successfully!", visible=True)
        }
        
    except Exception as e:
        current_pdb_info.update({"pdb_id": None, "pdb_path": None, "prepared_pdbqt": None, "docking_results": None, "prankweb_csv": None})
        return {
            info_box: gr.update(visible=False),
            structure_viewer: gr.update(value=""),
            download_file: gr.update(value=None),
            search_status: gr.update(value=f"❌ Error: {str(e)}", visible=True)
        }

# Function wrapper for ADMET to handle Gradio outputs
def process_admet():
    """
    Runs the ADMET pipeline and handles the output safely.
    Handles both success (tuple) and failure (None) from run_admet_prediction.
    """
    try:
        # Run prediction
        result = run_admet_prediction()
        
        # Handle Failure (e.g. no files found)
        if result is None:
            return {
                admet_status: gr.update(value="❌ Analysis Failed: No docking results found or SMILES generation failed. Please run docking first.", visible=True),
                admet_table: gr.update(visible=False),
                admet_download: gr.update(visible=False)
            }
            
        # Handle Success (unpack tuple)
        msg, df, csv_path = result
        
        return {
            admet_status: gr.update(value=f"✅ {msg}", visible=True),
            admet_table: gr.update(value=df, visible=True),
            admet_download: gr.update(value=csv_path, visible=True)
        }
        
    except Exception as e:
        # Unexpected System Error
        return {
            admet_status: gr.update(value=f"❌ System Error: {str(e)}", visible=True),
            admet_table: gr.update(visible=False),
            admet_download: gr.update(visible=False)
        }


# Create Gradio Interface with Tabs
with gr.Blocks(theme=gr.themes.Soft(), css="""
    .gradio-container {
        max-width: 1600px !important;
    }
    .main-header {
        text-align: center;
        padding: 40px 20px;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        border-radius: 20px;
        color: white;
        margin-bottom: 30px;
    }
    .main-header h1 {
        font-size: 42px;
        font-weight: 800;
        margin: 0 0 10px 0;
    }
    .main-header p {
        font-size: 18px;
        opacity: 0.95;
        margin: 0;
        font-weight: 500;
    }
    .nav-buttons {
        display: flex;
        justify-content: space-between;
        margin-top: 20px;
        padding-top: 20px;
        border-top: 2px solid #e5e7eb;
    }
""", title="Protein Structure Finder & Analyzer") as demo:
    
    gr.HTML("""
        <div class="main-header">
            <h1>🧬 Protein Structure Finder & Analyzer</h1>
            <p>Discover, visualize and analyze protein structures related to diseases</p>
        </div>
    """)
    
    with gr.Tabs() as tabs:
        # Tab 1: Structure Search
        with gr.Tab("🔍 Structure Search", id=0):
            with gr.Row():
                with gr.Column(scale=1):
                    disease_input = gr.Textbox(
                        label="🔍 Enter Disease or Condition",
                        placeholder="e.g., Alzheimer's Disease, diabetes, inflammation...",
                        lines=1
                    )
                    
                    search_btn = gr.Button("🚀 Search Structure", variant="primary", size="lg")
                    
                    info_box = gr.HTML(visible=False)
                    search_status = gr.Markdown(visible=False)
                    download_file = gr.File(label="💾 Download PDB File", visible=True)
                
                with gr.Column(scale=2):
                    structure_viewer = gr.HTML(label="🔬 3D Structure Viewer")
            
            with gr.Row(elem_classes="nav-buttons"):
                prev_btn_1 = gr.Button("← Previous", variant="secondary", visible=False)
                next_btn_1 = gr.Button("Next: Ramachandran Analysis →", variant="primary")
        
        # Tab 2: Ramachandran Plot
        with gr.Tab("📊 Ramachandran Analysis", id=1):
            gr.Markdown("### Ramachandran Plot Analysis")
            gr.Markdown("Analyze the backbone dihedral angles of your protein structure to assess its quality and validate the geometry.")
            
            ramplot_btn = gr.Button("🔬 Run Ramachandran Analysis", variant="secondary", size="lg")
            ramplot_status = gr.HTML(visible=False)
            
            with gr.Row():
                with gr.Column():
                    plot1 = gr.Image(label="Map Type 2D All", visible=False)
                with gr.Column():
                    plot2 = gr.Image(label="Map Type 3D All", visible=False)
            
            with gr.Row():
                with gr.Column():
                    plot3 = gr.Image(label="Std Map Type 2D General Gly", visible=False)
                with gr.Column():
                    plot4 = gr.Image(label="Std Map Type 3D General", visible=False)
            
            with gr.Row(elem_classes="nav-buttons"):
                prev_btn_2 = gr.Button("← Previous: Structure Search", variant="secondary")
                next_btn_2 = gr.Button("Next: Binding Site Prediction →", variant="primary")
        
        # Tab 3: PrankWeb Prediction
        with gr.Tab("🎯 Binding Site Prediction", id=2):
            gr.Markdown("### PrankWeb Binding Site Prediction")
            gr.Markdown("Predict potential ligand binding sites on your protein structure using PrankWeb.")
            
            prankweb_btn = gr.Button("🔮 Run PrankWeb Prediction", variant="secondary", size="lg")
            prankweb_status = gr.HTML(visible=False)
            prankweb_results = gr.Dataframe(label="Prediction Results", visible=False)
            
            with gr.Row(elem_classes="nav-buttons"):
                prev_btn_3 = gr.Button("← Previous: Ramachandran Analysis", variant="secondary")
                next_btn_3 = gr.Button("Next: Protein Preparation →", variant="primary")
        
        # Tab 4: Protein Preparation
        with gr.Tab("⚙️ Protein Preparation", id=3):
            gr.Markdown("### Protein Preparation for Docking (Meeko)")
            gr.Markdown("Prepare your protein structure for molecular docking by converting it to PDBQT format with proper charges and atom types.")
            
            prepare_btn = gr.Button("🔧 Prepare Protein with Meeko", variant="secondary", size="lg")
            prepare_status = gr.HTML(visible=False)
            
            with gr.Row():
                with gr.Column(scale=2):
                    prepared_viewer = gr.HTML(label="🔬 Prepared Structure Viewer")
                with gr.Column(scale=1):
                    prepared_download = gr.File(label="💾 Download PDBQT File", visible=True)
            
            with gr.Row(elem_classes="nav-buttons"):
                prev_btn_4 = gr.Button("← Previous: Binding Site Prediction", variant="secondary")
                next_btn_4 = gr.Button("Next: Molecular Docking →", variant="primary")
        
        # Tab 5: Molecular Docking
        with gr.Tab("🚀 Molecular Docking", id=4):
            gr.Markdown("### Molecular Docking (AutoDock Vina)")
            gr.Markdown("Perform molecular docking to predict how ligands bind to your protein target.")
            
            docking_btn = gr.Button("🚀 Run Molecular Docking", variant="secondary", size="lg")
            docking_status = gr.HTML(visible=False)
            docking_summary = gr.Dataframe(label="Docking Summary - Top 3 Poses per Ligand", visible=False)
            
            gr.Markdown("### 📊 View Docked Structures")
            
            with gr.Row():
                with gr.Column(scale=1):
                    pose_selector = gr.Dropdown(
                        label="Select Pose to View",
                        choices=[],
                        visible=False,
                        interactive=True
                    )
                    view_pose_btn = gr.Button("👁️ View Selected Pose", variant="primary", size="lg")
                
                with gr.Column(scale=2):
                    docked_viewer = gr.HTML(label="🔬 Docked Complex Viewer")
            
            with gr.Row(elem_classes="nav-buttons"):
                prev_btn_5 = gr.Button("← Previous: Protein Preparation", variant="secondary")
                next_btn_5 = gr.Button("Next: ADMET Analysis →", variant="primary")

        # Tab 6: ADMET Analysis
        with gr.Tab("🧪 ADMET Analysis", id=5):
            gr.Markdown("### 💊 Comprehensive ADMET Profiling & Filtering")
            gr.Markdown("""
            This module applies a **Multi-Stage Filtering Pipeline** to your docked ligands to ensure safety and developability:
            
            1.  **Primary Filters (Hard Stops):** Auto-rejects ligands with high toxicity risks (e.g., Ames Positive, High hERG/DILI).
            2.  **Developability Filters (Soft Constraints):** Flags ligands with suboptimal properties (e.g., MW > 550, SA Score > 6.0).
            3.  **Composite Scoring:** Calculates a final score (0-100) based on accumulated penalties.
            
            **Final Decision:** Ligands are classified as **ACCEPT**, **REVIEW**, or **REJECT** based on the rules.
            """)
            
            admet_btn = gr.Button("⚗️ Run ADMET Prediction", variant="secondary", size="lg")
            admet_status = gr.Markdown(visible=False)
            
            with gr.Row():
                admet_download = gr.File(label="📥 Download Full Report (CSV)", visible=False)
            
            admet_table = gr.Dataframe(
                label="ADMET Results Table", 
                visible=False,
                interactive=False,
                wrap=True
            )
            
            with gr.Row(elem_classes="nav-buttons"):
                prev_btn_6 = gr.Button("← Previous: Molecular Docking", variant="secondary")
                next_btn_6 = gr.Button("Back to Structure Search", variant="primary")
    
    # Navigation event handlers
    next_btn_1.click(lambda: gr.Tabs(selected=1), None, tabs)
    
    prev_btn_2.click(lambda: gr.Tabs(selected=0), None, tabs)
    next_btn_2.click(lambda: gr.Tabs(selected=2), None, tabs)
    
    prev_btn_3.click(lambda: gr.Tabs(selected=1), None, tabs)
    next_btn_3.click(lambda: gr.Tabs(selected=3), None, tabs)
    
    prev_btn_4.click(lambda: gr.Tabs(selected=2), None, tabs)
    next_btn_4.click(lambda: gr.Tabs(selected=4), None, tabs)
    
    prev_btn_5.click(lambda: gr.Tabs(selected=3), None, tabs)
    next_btn_5.click(lambda: gr.Tabs(selected=5), None, tabs)
    
    prev_btn_6.click(lambda: gr.Tabs(selected=4), None, tabs)
    next_btn_6.click(lambda: gr.Tabs(selected=0), None, tabs)
    
    # Main feature event handlers
    search_btn.click(
        fn=process_disease,
        inputs=[disease_input],
        outputs={info_box, structure_viewer, download_file, search_status}
    )
    
    ramplot_btn.click(
        fn=run_ramplot,
        inputs=[],
        outputs=[ramplot_status, plot1, plot2, plot3, plot4]
    )
    
    prankweb_btn.click(
        fn=run_prankweb_prediction,
        inputs=[],
        outputs=[prankweb_status, prankweb_results]
    )
    
    prepare_btn.click(
        fn=prepare_protein_meeko,
        inputs=[],
        outputs=[prepare_status, prepared_viewer, prepared_download]
    )
    
    docking_btn.click(
        fn=run_molecular_docking,
        inputs=[],
        outputs=[docking_status, docking_summary, pose_selector]
    )
    
    view_pose_btn.click(
        fn=display_docked_structure,
        inputs=[pose_selector],
        outputs=[docked_viewer]
    )

    # ADMET Handler
    admet_btn.click(
        fn=process_admet,
        inputs=[],
        outputs={admet_status, admet_table, admet_download}
    )

if __name__ == "__main__":
    demo.launch(share=False)