"""
Main Gradio interface for Protein Structure Finder & Analyzer
"""

import os
import gradio as gr
import pandas as pd

# Import modules
from config import current_pdb_info, PROTEINS_DIR
from ramachandran import run_ramplot
from prankweb import run_prankweb_prediction
from protein_prep import prepare_protein_meeko
from docking import run_molecular_docking  # Removed display_docked_structure as we define a new one here
from admet_analysis import run_admet_prediction
from utils import map_disease_to_protein, find_best_pdb_structure
from visualization import show_structure

def process_disease(user_input: str):
    """Main function to process disease/protein input."""
    
    if not user_input.strip():
        current_pdb_info.update({"pdb_id": None, "pdb_path": None})
        return {
            info_box: gr.update(visible=False),
            structure_viewer: gr.update(value=""),
            download_file: gr.update(value=None),
            search_status: gr.update(value="⚠️ Please enter a disease or protein name", visible=True)
        }
    
    protein_name = map_disease_to_protein(user_input)
    is_direct_protein = False
    
    if not protein_name:
        protein_name = user_input.strip()
        is_direct_protein = True
    
    result = find_best_pdb_structure(protein_name, max_check=100)
    
    if not result:
        current_pdb_info.update({"pdb_id": None, "pdb_path": None})
        error_msg = f"❌ No suitable PDB structure found for: {protein_name}"
        return {
            info_box: gr.update(visible=False),
            structure_viewer: gr.update(value=""),
            download_file: gr.update(value=None),
            search_status: gr.update(value=error_msg, visible=True)
        }
    
    pdb_id, pdb_path = result
    
    try:
        with open(pdb_path, 'r') as f:
            pdb_content = f.read()
        
        current_pdb_info.update({
            "pdb_id": pdb_id, 
            "pdb_path": pdb_path, 
            "prepared_pdbqt": None, 
            "docking_results": None, 
            "prankweb_csv": None
        })
        
        info_html = f"""
        <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 24px; border-radius: 16px; color: white;">
            <h3>Structure Loaded</h3>
            <p><strong>Input:</strong> {user_input}</p>
            <p><strong>Protein:</strong> {protein_name}</p>
            <p><strong>PDB ID:</strong> {pdb_id}</p>
            <p><strong>Species:</strong> Homo-sapiens</p>
        </div>
        """
        
        # UPDATED: Explicitly pass arguments to match new signature
        # show_structure(protein_text, ligand_text=None, pdb_id=..., protein_name=...)
        structure_html = show_structure(
            protein_text=pdb_content, 
            ligand_text=None, 
            pdb_id=pdb_id, 
            protein_name=protein_name
        )
        
        return {
            info_box: gr.update(value=info_html, visible=True),
            structure_viewer: gr.update(value=structure_html),
            download_file: gr.update(value=pdb_path),
            search_status: gr.update(value="✅ Structure loaded successfully!", visible=True)
        }
        
    except Exception as e:
        return {
            info_box: gr.update(visible=False),
            structure_viewer: gr.update(value=""),
            download_file: gr.update(value=None),
            search_status: gr.update(value=f"❌ Error: {str(e)}", visible=True)
        }

def visualize_docking_result(selection_value: str):
    """
    Visualizes the specific docked pose overlaid on the protein.
    Handles the 'filepath::pose_number' format correctly.
    """
    if not selection_value:
        return "⚠️ Please select a pose to view."
        
    try:
        # 1. Parse the selection string (Split path and pose number)
        if "::" not in selection_value:
            return f"❌ Invalid format. Expected 'path::pose', got: {selection_value}"
            
        ligand_path, pose_num_str = selection_value.split("::")
        target_pose_num = int(pose_num_str)
        
        # 2. Check if the actual files exist
        if not os.path.exists(ligand_path):
            return f"❌ Ligand file not found at: {ligand_path}"
            
        protein_path = current_pdb_info.get("pdb_path")
        if not protein_path or not os.path.exists(protein_path):
            return "❌ Protein structure not found. Please load a protein first."

        # 3. Read Protein Data
        with open(protein_path, 'r') as f:
            protein_text = f.read()
            
        # 4. Read Ligand Data & Extract Specific Pose
        # We need to find the block starting with "MODEL X" and ending with "ENDMDL"
        with open(ligand_path, 'r') as f:
            lines = f.readlines()
            
        model_lines = []
        in_model = False
        current_model = -1
        found_pose = False
        
        for line in lines:
            if line.startswith("MODEL"):
                try:
                    current_model = int(line.split()[1])
                except: pass
                if current_model == target_pose_num:
                    in_model = True
                    found_pose = True
            
            if in_model:
                model_lines.append(line)
                
            if line.startswith("ENDMDL") and in_model:
                in_model = False
                break 
        
        # Fallback: If we couldn't parse models (e.g. single pose file), use whole file
        ligand_text_pose = "".join(model_lines) if found_pose else "".join(lines)
        
        pdb_id = current_pdb_info.get("pdb_id", "Docking")
        
        # 5. Visualize
        return show_structure(
            protein_text=protein_text, 
            ligand_text=ligand_text_pose, 
            pdb_id=pdb_id, 
            protein_name=f"Docked Pose {target_pose_num}"
        )

    except Exception as e:
        return f"❌ Visualization Error: {str(e)}"

def process_admet():
    try:
        result = run_admet_prediction()
        if result is None:
            return {
                admet_status: gr.update(value="❌ Analysis Failed: Run docking first.", visible=True),
                admet_table: gr.update(visible=False),
                admet_download: gr.update(visible=False)
            }
        msg, df, csv_path = result
        return {
            admet_status: gr.update(value=f"✅ {msg}", visible=True),
            admet_table: gr.update(value=df, visible=True),
            admet_download: gr.update(value=csv_path, visible=True)
        }
    except Exception as e:
        return {
            admet_status: gr.update(value=f"❌ System Error: {str(e)}", visible=True),
            admet_table: gr.update(visible=False),
            admet_download: gr.update(visible=False)
        }

# UI Layout
with gr.Blocks(theme=gr.themes.Soft(), title="Protein Structure Finder & Analyzer") as demo:
    
    gr.HTML("<div class='main-header'><h1>🧬 Protein Structure Finder & Analyzer</h1></div>")
    
    with gr.Tabs() as tabs:
        # Tab 1: Search
        with gr.Tab("🔍 Structure Search", id=0):
            gr.Markdown("""
            ### Protein Identification
            * **Search:** Find PDB structures based on disease relevance or protein names.
            * **Validation:** View the experimental 3D structure and download the source PDB file.
            """)
            with gr.Row():
                with gr.Column(scale=1):
                    disease_input = gr.Textbox(label="Enter Disease/Gene/PDB ID", placeholder="e.g., Alzheimer's, Insulin")
                    search_btn = gr.Button("🚀 Search Best Structure", variant="primary")
                    info_box = gr.HTML(visible=False)
                    search_status = gr.Markdown(visible=False)
                    download_file = gr.File(label="Download PDB", visible=True)
                with gr.Column(scale=2):
                    structure_viewer = gr.HTML(label="3D Viewer")
            
            with gr.Row():
                next_btn_1 = gr.Button("Next: Ramachandran Analysis →", variant="primary")
        
        # Tab 2: Ramachandran
        with gr.Tab("📊 Ramachandran Analysis", id=1):
            gr.Markdown("""
            ### Structural Quality Assessment
            * **Geometry Check:** Evaluates backbone dihedral angles ($\phi$ and $\psi$) to ensure structural validity.
            * **Automated Repair:** If significant residues are missing, the system uses homology modeling (Swiss-Model) to refine the structure.
            """)
            ramplot_btn = gr.Button("🔬 Run Ramachandran Analysis", variant="secondary")
            
            ramplot_status = gr.HTML(visible=False)
            
            with gr.Row():
                plot1 = gr.Image(label="Map Type 2D", visible=False)
                plot2 = gr.Image(label="Map Type 3D", visible=False)
            with gr.Row():
                plot3 = gr.Image(label="Std Map 2D", visible=False)
                plot4 = gr.Image(label="Std Map 3D", visible=False)
            
            with gr.Row():
                prev_btn_2 = gr.Button("← Previous", variant="secondary")
                next_btn_2 = gr.Button("Next: Protein Preparation →", variant="primary")
        
        # Tab 3: Protein Preparation
        with gr.Tab("⚙️ Protein Preparation", id=2):
            gr.Markdown("""
            ### File Preparation for Simulation
            * **PDBQT Conversion:** Prepares the protein for docking engines like AutoDock Vina.
            * **Refinement:** Removes water molecules, adds polar hydrogens, and computes partial charges.
            """)
            prepare_btn = gr.Button("🔧 Prepare Protein", variant="secondary")
            prepare_status = gr.HTML(visible=False)
            with gr.Row():
                prepared_viewer = gr.HTML(label="Prepared Viewer")
                prepared_download = gr.File(label="Download PDBQT")
            with gr.Row():
                prev_btn_3 = gr.Button("← Previous", variant="secondary")
                next_btn_3 = gr.Button("Next: Binding Site Prediction →", variant="primary")

        # Tab 4: PrankWeb
        with gr.Tab("🎯 Binding Site Prediction", id=3):
            gr.Markdown("""
            ### Active Site Prediction
            * **Pocket Identification:** Uses machine learning to detect potential binding cavities on the protein surface.
            * **Guidance:** Helps define the search space for molecular docking.
            """)
            prankweb_btn = gr.Button("🔮 Run PrankWeb", variant="secondary")
            prankweb_status = gr.HTML(visible=False)
            prankweb_results = gr.Dataframe(label="Results", visible=False)
            with gr.Row():
                prev_btn_4 = gr.Button("← Previous", variant="secondary")
                next_btn_4 = gr.Button("Next: Docking →", variant="primary")

        # Tab 5: Docking
        with gr.Tab("🚀 Molecular Docking", id=4):
            gr.Markdown("""
            ### Molecular Docking
            * **Interaction Study:** Simulates how a ligand binds to the target protein.
            * **Scoring:** Calculates binding affinity (kcal/mol) for multiple docking poses.
            """)
            docking_btn = gr.Button("Run Docking", variant="secondary")
            docking_status = gr.HTML(visible=False)
            docking_summary = gr.Dataframe(visible=False)
            
            # The dropdown will be populated by run_molecular_docking with file paths
            pose_selector = gr.Dropdown(label="Select Pose to View", visible=False)
            
            view_pose_btn = gr.Button("View Pose", variant="primary")
            docked_viewer = gr.HTML(label="Docked Interaction Viewer")
            
            with gr.Row():
                prev_btn_5 = gr.Button("← Previous", variant="secondary")
                next_btn_5 = gr.Button("Next: ADMET →", variant="primary")

        # Tab 6: ADMET
        with gr.Tab("🧪 ADMET Analysis", id=5):
            gr.Markdown("""
            ### Drug-likeness & Safety
            * **Pharmacokinetics:** Evaluates Absorption, Distribution, Metabolism, and Excretion.
            * **Toxicity:** Assesses the potential risk and safety profile of the molecule.
            """)
            admet_btn = gr.Button("Run ADMET", variant="secondary")
            admet_status = gr.Markdown(visible=False)
            admet_download = gr.File(visible=False)
            admet_table = gr.Dataframe(visible=False)
            with gr.Row():
                prev_btn_6 = gr.Button("← Previous", variant="secondary")
                next_btn_6 = gr.Button("Back to Start", variant="primary")

    # Events
    next_btn_1.click(lambda: gr.Tabs(selected=1), None, tabs)
    
    prev_btn_2.click(lambda: gr.Tabs(selected=0), None, tabs)
    next_btn_2.click(lambda: gr.Tabs(selected=2), None, tabs)
    
    # Prep Tab Events
    prev_btn_3.click(lambda: gr.Tabs(selected=1), None, tabs)
    next_btn_3.click(lambda: gr.Tabs(selected=3), None, tabs)
    
    # PrankWeb Tab Events
    prev_btn_4.click(lambda: gr.Tabs(selected=2), None, tabs)
    next_btn_4.click(lambda: gr.Tabs(selected=4), None, tabs)
    
    prev_btn_5.click(lambda: gr.Tabs(selected=3), None, tabs)
    next_btn_5.click(lambda: gr.Tabs(selected=5), None, tabs)
    
    prev_btn_6.click(lambda: gr.Tabs(selected=4), None, tabs)
    next_btn_6.click(lambda: gr.Tabs(selected=0), None, tabs)

    search_btn.click(process_disease, inputs=[disease_input], outputs={info_box, structure_viewer, download_file, search_status})
    ramplot_btn.click(fn=run_ramplot, inputs=[], outputs=[ramplot_status, plot1, plot2, plot3, plot4])
    prankweb_btn.click(fn=run_prankweb_prediction, inputs=[], outputs=[prankweb_status, prankweb_results])
    prepare_btn.click(fn=prepare_protein_meeko, inputs=[], outputs=[prepare_status, prepared_viewer, prepared_download])
    docking_btn.click(fn=run_molecular_docking, inputs=[], outputs=[docking_status, docking_summary, pose_selector])
    
    # UPDATED: Use the new local visualization wrapper
    view_pose_btn.click(fn=visualize_docking_result, inputs=[pose_selector], outputs=[docked_viewer])
    
    admet_btn.click(fn=process_admet, inputs=[], outputs={admet_status, admet_table, admet_download})

if __name__ == "__main__":
    demo.launch(share=True)