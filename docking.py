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
    
    # Check if Vina executable exists
    if not os.path.exists(VINA_EXE) and not shutil.which(VINA_EXE):
         # Fallback check if it's in path, otherwise warn
         # For this specific script, we assume it's in the root or defined path
         pass 

    # Show processing message
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
        
        # First, try to use the stored CSV path
        if current_pdb_info.get("prankweb_csv") and os.path.exists(current_pdb_info["prankweb_csv"]):
            csv_file = current_pdb_info["prankweb_csv"]
        else:
            # Search for the CSV file in all subdirectories
            for root, dirs, files in os.walk(output_dir):
                for file in files:
                    if file.endswith("_predictions.csv"):
                        csv_file = os.path.join(root, file)
                        break
                if csv_file:
                    break
        
        # Check if PrankWeb results exist
        if not csv_file or not os.path.exists(csv_file):
            yield (
                gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ PrankWeb results not found. Please run PrankWeb prediction first.</div>", visible=True),
                gr.update(value=None, visible=False),
                gr.update(choices=[], visible=False)
            )
            return
        
        # Check if ligand folder exists
        if not os.path.exists(ligand_folder):
            yield (
                gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Ligand folder 'pdbqt' not found. Please ensure ligands are prepared.</div>", visible=True),
                gr.update(value=None, visible=False),
                gr.update(choices=[], visible=False)
            )
            return
        
        # Load pocket table
        df = pd.read_csv(csv_file)
        
        # Get all ligand PDBQT files
        ligand_files = glob.glob(os.path.join(ligand_folder, "*.pdbqt"))
        
        if not ligand_files:
            yield (
                gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ No ligand files found in 'pdbqt' folder.</div>", visible=True),
                gr.update(value=None, visible=False),
                gr.update(choices=[], visible=False)
            )
            return
        
        # Ensure output folders exist
        output_dir_pdbqt = os.path.join(DOCKING_RESULTS_DIR, "pdbqt")
        output_dir_pdb = os.path.join(DOCKING_RESULTS_DIR, "pdb")
        os.makedirs(output_dir_pdbqt, exist_ok=True)
        os.makedirs(output_dir_pdb, exist_ok=True)
        
        # Initialize summary data
        summary_data = []
        
        # Helper functions
        def convert_pdbqt_to_pdb(pdbqt_file, pdb_file):
            try:
                subprocess.run(['obabel', pdbqt_file, '-O', pdb_file, '-h'], 
                             check=True, capture_output=True)
                return True
            except:
                return False
        
        def save_protein_ligand_complex(protein_pdbqt, ligand_poses_pdbqt, output_pdb, ligand_name, pocket_name):
            try:
                protein_pdb = f"temp_protein_{ligand_name}_{pocket_name}.pdb"
                subprocess.run(['obabel', protein_pdbqt, '-O', protein_pdb, '-h'], 
                             check=True, capture_output=True)
                
                ligand_pdb = f"temp_ligand_{ligand_name}_{pocket_name}.pdb"
                subprocess.run(['obabel', ligand_poses_pdbqt, '-O', ligand_pdb, '-h'], 
                             check=True, capture_output=True)
                
                with open(output_pdb, 'w') as outfile:
                    with open(protein_pdb, 'r') as prot:
                        outfile.write(prot.read())
                    with open(ligand_pdb, 'r') as lig:
                        outfile.write(lig.read())
                
                os.remove(protein_pdb)
                os.remove(ligand_pdb)
                return True
            except:
                return False
        
        # Iterate over each ligand file
        for ligand_pdbqt in ligand_files:
            ligand_name = os.path.splitext(os.path.basename(ligand_pdbqt))[0]
            ligand_best_poses = []
            
            # Iterate over each pocket for this ligand
            for index, row in df.iterrows():
                pocket_name = row['name     '].strip()
                center_x = float(row['   center_x'])
                center_y = float(row['   center_y'])
                center_z = float(row['   center_z'])
                
                # Define output file for Vina
                output_pdbqt_file = os.path.join(output_dir_pdbqt, f"{ligand_name}_{pocket_name}_docked_poses.pdbqt")

                # Construct Vina command
                cmd = [
                    VINA_EXE,
                    "--receptor", protein_pdbqt,
                    "--ligand", ligand_pdbqt,
                    "--center_x", str(center_x),
                    "--center_y", str(center_y),
                    "--center_z", str(center_z),
                    "--size_x", "25",
                    "--size_y", "25",
                    "--size_z", "25",
                    "--exhaustiveness", "8",
                    "--num_modes", "10",
                    "--out", output_pdbqt_file
                ]

                # Run Vina via subprocess
                try:
                    result = subprocess.run(
                        cmd, capture_output=True, text=True, check=True, shell=True
                    )
                    
                    # Parse binding energies from stdout
                    # Vina output format: 
                    # mode |   affinity | dist from best mode
                    #      | (kcal/mol) | rmsd l.b.| rmsd u.b.
                    # -----+------------+----------+----------
                    #    1       -7.5      0.000      0.000
                    
                    for line in result.stdout.splitlines():
                        # Regex to find lines starting with a number followed by float (energy)
                        if re.match(r'^\s*\d+', line):
                            # Extract parts: mode, affinity, ...
                            parts = line.split()
                            if len(parts) >= 2:
                                try:
                                    mode_num = int(parts[0])
                                    affinity = float(parts[1])
                                    
                                    # Filter poses with binding energy <= -7.0 kcal/mol
                                    if affinity <= -7.0:
                                        ligand_best_poses.append({
                                            'ligand': ligand_name,
                                            'pocket': pocket_name,
                                            'pose_number': mode_num,
                                            'binding_energy': affinity,
                                            'center_x': center_x,
                                            'center_y': center_y,
                                            'center_z': center_z
                                        })
                                except ValueError:
                                    continue

                    # Proceed to conversion only if output file was created
                    if os.path.exists(output_pdbqt_file):
                        # Save ligand-only PDB
                        pdb_ligand_only = os.path.join(output_dir_pdb, f"{ligand_name}_{pocket_name}_ligand_poses.pdb")
                        convert_pdbqt_to_pdb(output_pdbqt_file, pdb_ligand_only)
                        
                        # Save protein + ligand complex PDB
                        pdb_complex = os.path.join(output_dir_pdb, f"{ligand_name}_{pocket_name}_complex.pdb")
                        save_protein_ligand_complex(protein_pdbqt, output_pdbqt_file, pdb_complex, ligand_name, pocket_name)
                        
                except subprocess.CalledProcessError as e:
                    print(f"Docking failed for {ligand_name} in {pocket_name}: {e}")
                    continue
            
            # Select top 3 poses for this ligand (across all pockets)
            if ligand_best_poses:
                ligand_best_poses.sort(key=lambda x: x['binding_energy'])
                top_3 = ligand_best_poses[:3]
                summary_data.extend(top_3)
        
        # Create summary table
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            summary_df = summary_df[['ligand', 'pocket', 'pose_number', 'binding_energy', 
                                    'center_x', 'center_y', 'center_z']]
            summary_df['binding_energy'] = summary_df['binding_energy'].round(3)
            summary_df['center_x'] = summary_df['center_x'].round(2)
            summary_df['center_y'] = summary_df['center_y'].round(2)
            summary_df['center_z'] = summary_df['center_z'].round(2)
            
            # Save summary table
            summary_file = os.path.join(DOCKING_RESULTS_DIR, "docking_summary.csv")
            summary_df.to_csv(summary_file, index=False)
            
            # Store in global variable
            current_pdb_info["docking_results"] = summary_df
            
            # Create dropdown choices
            choices = []
            for idx, row in summary_df.iterrows():
                label = f"{row['ligand']} - {row['pocket']} (Pose {row['pose_number']}) | {row['binding_energy']:.2f} kcal/mol"
                choices.append(label)
            
            success_msg = "<div style='padding: 20px; background: #d4edda; border-radius: 8px; color: #155724;'>"
            success_msg += f"✅ Docking completed!<br>"
            success_msg += f"<small>Processed {len(ligand_files)} ligands</small><br>"
            success_msg += f"<small>Found {len(summary_data)} good poses (≤ -7.0 kcal/mol)</small><br>"
            success_msg += f"<small>Best energy: {summary_df['binding_energy'].min():.2f} kcal/mol</small>"
            success_msg += "</div>"
            
            yield (
                gr.update(value=success_msg, visible=True),
                gr.update(value=summary_df, visible=True),
                gr.update(choices=choices, visible=True, value=choices[0] if choices else None)
            )
        else:
            yield (
                gr.update(value="<div style='padding: 20px; background: #f8d7da; border-radius: 8px; color: #721c24;'>⚠️ No poses with binding energy ≤ -7.0 kcal/mol found</div>", visible=True),
                gr.update(value=None, visible=False),
                gr.update(choices=[], visible=False)
            )
    
    except Exception as e:
        yield (
            gr.update(value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Error: {str(e)}</div>", visible=True),
            gr.update(value=None, visible=False),
            gr.update(choices=[], visible=False)
        )


def display_docked_structure(selected_pose):
    """Display the selected docked structure in 3D."""
    if not selected_pose:
        return gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Please select a pose to display</div>")
    
    if "docking_results" not in current_pdb_info:
        return gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ No docking results available</div>")
    
    try:
        # Parse the selected pose label
        # Format: "ligand - pocket (Pose N) | energy kcal/mol"
        parts = selected_pose.split(" - ")
        ligand_name = parts[0].strip()
        
        pocket_part = parts[1].split(" (Pose ")[0].strip()
        pose_num_part = parts[1].split(" (Pose ")[1].split(")")[0].strip()
        
        # Find the complex PDB file
        complex_file = os.path.join(DOCKING_RESULTS_DIR, "pdb", f"{ligand_name}_{pocket_part}_complex.pdb")
        
        if not os.path.exists(complex_file):
            return gr.update(value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Complex file not found: {complex_file}</div>")
        
        # Read PDB content
        with open(complex_file, 'r') as f:
            pdb_content = f.read()
        
        # Create 3D visualization
        structure_html = show_structure(pdb_content, f"{ligand_name}_{pocket_part}", f"Docked Complex - {selected_pose}")
        
        return gr.update(value=structure_html)
        
    except Exception as e:
        return gr.update(value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Error displaying structure: {str(e)}</div>")