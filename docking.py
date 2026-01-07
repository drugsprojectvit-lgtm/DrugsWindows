"""
Molecular docking module using AutoDock Vina Executable
"""

import os
import glob
import subprocess
import re
import pandas as pd
import gradio as gr
from config import current_pdb_info, DOCKING_RESULTS_DIR, LIGAND_DIR, PRANKWEB_OUTPUT_DIR
from visualization import show_structure

# Path to your Vina executable
VINA_EXE = "vina.exe" 

def run_molecular_docking():
    """Run molecular docking using AutoDock Vina (executable via subprocess)."""
    
    if not current_pdb_info.get("prepared_pdbqt"):
        return (
            gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ No prepared protein found. Please prepare protein first.</div>", visible=True),
            gr.update(value=None, visible=False),
            gr.update(choices=[], visible=False)
        )
    
    yield (
        gr.update(value="<div style='padding: 20px; background: #fff3cd; border-radius: 8px; color: #856404;'>🔬 Running molecular docking (this may take several minutes)...</div>", visible=True),
        gr.update(value=None, visible=False),
        gr.update(choices=[], visible=False)
    )
    
    try:
        # Input files and directories
        protein_pdbqt = current_pdb_info["prepared_pdbqt"]
        ligand_folder = LIGAND_DIR
        
        # Dynamically find the CSV file from PrankWeb results
        output_dir = PRANKWEB_OUTPUT_DIR
        csv_file = None
        
        if current_pdb_info.get("prankweb_csv") and os.path.exists(current_pdb_info["prankweb_csv"]):
            csv_file = current_pdb_info["prankweb_csv"]
        else:
            for root, dirs, files in os.walk(output_dir):
                for file in files:
                    if file.endswith("_predictions.csv"):
                        csv_file = os.path.join(root, file)
                        break
                if csv_file: break
        
        # Validation Checks
        if not csv_file or not os.path.exists(csv_file):
            yield (gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ PrankWeb results not found.</div>", visible=True), None, gr.update(choices=[]))
            return
        
        if not os.path.exists(ligand_folder):
            yield (gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Ligand folder not found.</div>", visible=True), None, gr.update(choices=[]))
            return
        
        # --- FIX: ROBUST CSV LOADING ---
        df = pd.read_csv(csv_file)
        
        # Strip whitespace from column names to fix 'name    ' vs 'name' issues
        df.columns = df.columns.str.strip()
        
        # Verify required columns exist
        required_cols = ['name', 'center_x', 'center_y', 'center_z']
        missing = [col for col in required_cols if col not in df.columns]
        if missing:
            yield (gr.update(value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ CSV format error. Missing columns: {missing}</div>", visible=True), None, gr.update(choices=[]))
            return
        # -------------------------------

        ligand_files = glob.glob(os.path.join(ligand_folder, "*.pdbqt"))
        
        if not ligand_files:
            yield (gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ No ligand files found.</div>", visible=True), None, gr.update(choices=[]))
            return
        
        # Directories
        output_dir_pdbqt = os.path.join(DOCKING_RESULTS_DIR, "pdbqt")
        output_dir_pdb = os.path.join(DOCKING_RESULTS_DIR, "pdb")
        os.makedirs(output_dir_pdbqt, exist_ok=True)
        os.makedirs(output_dir_pdb, exist_ok=True)
        
        summary_data = []
        
        # --- Processing Loop ---
        for ligand_pdbqt in ligand_files:
            ligand_name = os.path.splitext(os.path.basename(ligand_pdbqt))[0]
            ligand_best_poses = []
            
            for index, row in df.iterrows():
                # Use CLEAN column names
                pocket_name = str(row['name']).strip()
                cx, cy, cz = float(row['center_x']), float(row['center_y']), float(row['center_z'])
                
                output_pdbqt_file = os.path.join(output_dir_pdbqt, f"{ligand_name}_{pocket_name}_docked_poses.pdbqt")

                # Vina Command
                cmd = [
                    VINA_EXE, "--receptor", protein_pdbqt, "--ligand", ligand_pdbqt,
                    "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
                    "--size_x", "25", "--size_y", "25", "--size_z", "25",
                    "--exhaustiveness", "8", "--num_modes", "10", "--out", output_pdbqt_file
                ]

                try:
                    result = subprocess.run(cmd, capture_output=True, text=True, check=True, shell=True)
                    
                    # Parse Energies
                    for line in result.stdout.splitlines():
                        if re.match(r'^\s*\d+', line):
                            parts = line.split()
                            if len(parts) >= 2:
                                try:
                                    mode_num = int(parts[0])
                                    affinity = float(parts[1])
                                    
                                    # Filter good poses
                                    if affinity <= -7.0:
                                        ligand_best_poses.append({
                                            'ligand': ligand_name,
                                            'pocket': pocket_name,
                                            'pose_number': mode_num,
                                            'binding_energy': affinity,
                                            'center_x': cx, 'center_y': cy, 'center_z': cz,
                                            # Store path to ligand-only PDB
                                            'pdb_file': os.path.join(output_dir_pdb, f"{ligand_name}_{pocket_name}_ligand_poses.pdb")
                                        })
                                except ValueError: continue

                    # 1. Convert PDBQT -> PDB (Ligand Only) via Obabel
                    if os.path.exists(output_pdbqt_file):
                        pdb_ligand_only = os.path.join(output_dir_pdb, f"{ligand_name}_{pocket_name}_ligand_poses.pdb")
                        # We use -h to add hydrogens back for better visualization
                        subprocess.run(['obabel', output_pdbqt_file, '-O', pdb_ligand_only, '-h'], check=False, capture_output=True)

                        # --- NEW: Generate Combined Complex (Protein + Ligand HETATM) ---
                        try:
                            complex_file = os.path.join(output_dir_pdb, f"{ligand_name}_{pocket_name}_complex.pdb")
                            
                            # A. Get Receptor Lines (from original PDB to preserve headers/chains if available)
                            receptor_path = current_pdb_info.get("pdb_path")
                            receptor_lines = []
                            if receptor_path and os.path.exists(receptor_path):
                                with open(receptor_path, 'r') as rf:
                                    # Keep ATOM, HETATM (cofactors), and TER records
                                    receptor_lines = [l for l in rf if l.startswith(("ATOM", "HETATM", "TER"))]
                            
                            # B. Get Ligand Lines (Extract BEST POSE/MODEL 1 from Vina output)
                            ligand_lines = []
                            with open(output_pdbqt_file, 'r') as lf:
                                in_model_1 = False
                                for line in lf:
                                    if line.startswith("MODEL 1"):
                                        in_model_1 = True
                                        continue
                                    if line.startswith("ENDMDL"):
                                        # Only grab the best pose (Model 1) for the complex file
                                        break 
                                    
                                    if in_model_1 and line.startswith("ATOM"):
                                        # CHANGE: Rename ATOM to HETATM for differentiation
                                        new_line = "HETATM" + line[6:]
                                        ligand_lines.append(new_line)

                            # C. Write Combined File
                            with open(complex_file, 'w') as cf:
                                cf.writelines(receptor_lines)
                                # Ensure separation
                                if receptor_lines and not receptor_lines[-1].strip() == "TER":
                                    cf.write("TER\n")
                                cf.writelines(ligand_lines)
                                cf.write("END\n")

                        except Exception as complex_err:
                            print(f"Warning: Could not generate complex PDB for {ligand_name}: {complex_err}")
                        # ----------------------------------------------------------------

                except subprocess.CalledProcessError as e:
                    print(f"Docking failed for {ligand_name}: {e}")
                    continue
            
            # Keep top 3 poses per ligand
            if ligand_best_poses:
                ligand_best_poses.sort(key=lambda x: x['binding_energy'])
                summary_data.extend(ligand_best_poses[:3])
        
        # --- Summary & Output ---
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            
            choices = []
            for idx, row in summary_df.iterrows():
                label = f"{row['ligand']} - {row['pocket']} (Pose {row['pose_number']}) | {row['binding_energy']:.2f} kcal/mol"
                value = f"{row['pdb_file']}::{row['pose_number']}"
                choices.append((label, value))
            
            success_msg = f"""<div style='padding: 20px; background: #d4edda; border-radius: 8px; color: #155724;'>
                ✅ Docking completed! Found {len(summary_data)} poses ≤ -7.0 kcal/mol.
            </div>"""
            
            yield (
                gr.update(value=success_msg, visible=True),
                gr.update(value=summary_df, visible=True),
                gr.update(choices=choices, visible=True, value=choices[0][1] if choices else None)
            )
        else:
            yield (
                gr.update(value="<div style='padding: 20px; background: #f8d7da; border-radius: 8px; color: #721c24;'>⚠️ No good poses found (≤ -7.0 kcal/mol)</div>", visible=True),
                gr.update(value=None, visible=False),
                gr.update(choices=[], visible=False)
            )
    
    except Exception as e:
        yield (gr.update(value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Error: {str(e)}</div>", visible=True), None, gr.update(choices=[]))


def display_docked_structure(selection_value):
    """
    Displays the docking result by overlaying the Protein and the specific Ligand Pose.
    """
    if not selection_value:
        return "<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Please select a pose.</div>"
    
    try:
        if "::" not in selection_value:
             return f"❌ Invalid selection format: {selection_value}"
             
        ligand_path, pose_num_str = selection_value.split("::")
        target_pose_num = int(pose_num_str)
        
        # Get Protein Structure
        protein_path = current_pdb_info.get("pdb_path")
        if not protein_path or not os.path.exists(protein_path):
             return "❌ Protein structure not found in memory."
             
        with open(protein_path, 'r') as f:
            protein_text = f.read()

        # Get Ligand Structure (Specific Pose)
        if not os.path.exists(ligand_path):
             return f"❌ Ligand file not found: {ligand_path}"

        # Extract only the specific model
        ligand_text_pose = ""
        with open(ligand_path, 'r') as f:
            lines = f.readlines()
            
        in_model = False
        current_model_num = -1
        found_pose = False
        model_lines = []
        
        for line in lines:
            if line.startswith("MODEL"):
                try:
                    current_model_num = int(line.split()[1])
                except: pass
                if current_model_num == target_pose_num:
                    in_model = True
                    found_pose = True
            
            if in_model:
                model_lines.append(line)
            
            if line.startswith("ENDMDL"):
                if in_model:
                    in_model = False
                    break 
        
        if found_pose:
            ligand_text_pose = "".join(model_lines)
        else:
            # Fallback if no MODEL tags (single pose file)
            ligand_text_pose = "".join(lines)

        return show_structure(
            protein_text=protein_text,
            ligand_text=ligand_text_pose,
            pdb_id="Docking Result",
            protein_name=f"Pose {target_pose_num}"
        )
        
    except Exception as e:
        return f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Visualization Error: {str(e)}</div>"