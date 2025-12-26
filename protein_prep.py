"""
Protein preparation module using Meeko
"""

import os
import subprocess
import tempfile
import gradio as gr
from config import current_pdb_info, PREPARED_PROTEIN_DIR
from visualization import show_structure


def prepare_protein_meeko():
    """Prepare protein using Meeko for docking."""
    if not current_pdb_info["pdb_id"] or not current_pdb_info["pdb_path"]:
        return (
            gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ No structure loaded. Please search for a disease first.</div>", visible=True),
            gr.update(value=""),
            gr.update(value=None)
        )
    
    # Show processing message
    yield (
        gr.update(value="<div style='padding: 20px; background: #fff3cd; border-radius: 8px; color: #856404;'>⚙️ Preparing protein with Meeko...</div>", visible=True),
        gr.update(value=""),
        gr.update(value=None)
    )
    
    pdb_path = current_pdb_info["pdb_path"]
    pdb_id = current_pdb_info["pdb_id"]
    
    output_dir = PREPARED_PROTEIN_DIR
    os.makedirs(output_dir, exist_ok=True)
    
    output_base = os.path.join(output_dir, "prepared_protein")
    
    cmd = [
        'mk_prepare_receptor.py',
        '-i', pdb_path,
        '-o', output_base,
        '-p',
        '--charge_model', 'gasteiger',
        '--default_altloc', 'A'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # The output file will be output_base.pdbqt
        pdbqt_path = f"{output_base}.pdbqt"
        
        if not os.path.exists(pdbqt_path):
            yield (
                gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ PDBQT file not generated</div>", visible=True),
                gr.update(value=""),
                gr.update(value=None)
            )
            return
        
        # Store prepared protein path globally
        current_pdb_info["prepared_pdbqt"] = pdbqt_path
        
        # Read PDBQT content
        with open(pdbqt_path, 'r') as f:
            pdbqt_content = f.read()
        
        # Create 3D visualization (PDBQT format is similar to PDB)
        protein_name = f"Prepared Protein ({pdb_id})"
        structure_html = show_structure(pdbqt_content, pdb_id, protein_name)
        
        # Create download file
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.pdbqt', delete=False)
        temp_file.write(pdbqt_content)
        temp_file.close()
        
        success_msg = "<div style='padding: 20px; background: #d4edda; border-radius: 8px; color: #155724;'>"
        success_msg += "✅ Protein preparation completed!<br>"
        success_msg += f"<small>Output: {pdbqt_path}</small>"
        if result.stdout:
            success_msg += f"<br><small>{result.stdout}</small>"
        success_msg += "</div>"
        
        yield (
            gr.update(value=success_msg, visible=True),
            gr.update(value=structure_html),
            gr.update(value=temp_file.name)
        )
        
    except subprocess.CalledProcessError as e:
        error_msg = f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>"
        error_msg += f"⚠️ Preparation failed:<br><small>{e.stderr if e.stderr else str(e)}</small></div>"
        yield (
            gr.update(value=error_msg, visible=True),
            gr.update(value=""),
            gr.update(value=None)
        )
    except Exception as e:
        yield (
            gr.update(value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Error: {str(e)}</div>", visible=True),
            gr.update(value=""),
            gr.update(value=None)
        )

