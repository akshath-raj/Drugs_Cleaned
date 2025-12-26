"""
Ramachandran plot analysis module
"""

import os
import subprocess
import gradio as gr
from config import current_pdb_info, RAMPLOT_OUTPUT_DIR, PROTEINS_DIR


def run_ramplot():
    """Run RAMPlot analysis on the current PDB file."""
    if not current_pdb_info["pdb_id"] or not current_pdb_info["pdb_path"]:
        return (
            gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ No structure loaded. Please search for a disease first.</div>", visible=True),
            gr.update(value=None, visible=False),
            gr.update(value=None, visible=False),
            gr.update(value=None, visible=False),
            gr.update(value=None, visible=False)
        )
    
    # Show processing message
    yield (
        gr.update(value="<div style='padding: 20px; background: #fff3cd; border-radius: 8px; color: #856404;'>🔬 Processing Ramachandran plot analysis...</div>", visible=True),
        gr.update(value=None, visible=False),
        gr.update(value=None, visible=False),
        gr.update(value=None, visible=False),
        gr.update(value=None, visible=False)
    )
    
    pdb_id = current_pdb_info["pdb_id"]
    input_folder = PROTEINS_DIR
    output_folder = RAMPLOT_OUTPUT_DIR
    
    os.makedirs(input_folder, exist_ok=True)
    os.makedirs(output_folder, exist_ok=True)
    
    cmd = [
        "ramplot", "pdb",
        "-i", input_folder,
        "-o", output_folder,
        "-m", "0",
        "-r", "600",
        "-p", "png"
    ]
    
    try:
        result = subprocess.run(cmd, check=True, text=True, capture_output=True)
        
        # Look for the generated plot files
        plot_files = {
            'map2d': os.path.join(output_folder, "Plots", "MapType2DAll.png"),
            'map3d': os.path.join(output_folder, "Plots", "MapType3DAll.png"),
            'std2d': os.path.join(output_folder, "Plots", "StdMapType2DGeneralGly.png"),
            'std3d': os.path.join(output_folder, "Plots", "StdMapType3DGeneral.png")
        }
        
        yield (
            gr.update(value="<div style='padding: 20px; background: #d4edda; border-radius: 8px; color: #155724;'>✅ Ramachandran plot analysis completed!</div>", visible=True),
            gr.update(value=plot_files['map2d'], visible=True),
            gr.update(value=plot_files['map3d'], visible=True),
            gr.update(value=plot_files['std2d'], visible=True),
            gr.update(value=plot_files['std3d'], visible=True)
        )
        
    except subprocess.CalledProcessError as e:
        yield (
            gr.update(value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>⚠️ Analysis failed: {e.stderr}</div>", visible=True),
            gr.update(value=None, visible=False),
            gr.update(value=None, visible=False),
            gr.update(value=None, visible=False),
            gr.update(value=None, visible=False)
        )
    except Exception as e:
        yield (
            gr.update(value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Error: {str(e)}</div>", visible=True),
            gr.update(value=None, visible=False),
            gr.update(value=None, visible=False),
            gr.update(value=None, visible=False),
            gr.update(value=None, visible=False)
        )
