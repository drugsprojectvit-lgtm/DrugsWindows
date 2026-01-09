"""
PrankWeb (Local P2Rank) binding site prediction module
"""

import os
import subprocess
import glob
import pandas as pd
import gradio as gr
from config import current_pdb_info, PRANKWEB_OUTPUT_DIR, PREPARED_PROTEIN_DIR

# --- CONFIGURATION ---
# Path to your P2Rank folder (relative or absolute)
P2RANK_HOME = "p2rank_2.5.1" 
# On Windows use 'prank.bat', on Linux/Mac use 'prank'
P2RANK_CMD = "prank.bat" if os.name == 'nt' else "prank"

def convert_pdbqt_to_pdb(pdbqt_file, pdb_file):
    """
    Converts PDBQT to PDB using OpenBabel (obabel).
    """
    try:
        # Using obabel to convert. -h adds hydrogens
        subprocess.run(['obabel', pdbqt_file, '-O', pdb_file], 
                       check=True, capture_output=True, shell=True)
        return True
    except Exception as e:
        print(f"OpenBabel conversion failed: {e}")
        return False

def run_prankweb_prediction():
    """Run Binding Site prediction using LOCAL P2Rank."""
    
    # 1. Validate Input (Must come from Protein Prep step)
    prepared_pdbqt_path = current_pdb_info.get("prepared_pdbqt")
    
    # Fallback checks
    if not prepared_pdbqt_path:
         potential_path = os.path.join(PREPARED_PROTEIN_DIR, "prepared_protein.pdbqt")
         if os.path.exists(potential_path):
             prepared_pdbqt_path = potential_path

    if not prepared_pdbqt_path or not os.path.exists(prepared_pdbqt_path):
        return (
            gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ No prepared protein found. Please run <b>Protein Preparation</b> first.</div>", visible=True),
            gr.update(value=None, visible=False)
        )
    
    # Check if P2Rank exists
    if not os.path.exists(P2RANK_HOME):
         return (
            gr.update(value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ P2Rank folder not found at: {P2RANK_HOME}. Please check the path.</div>", visible=True),
            gr.update(value=None, visible=False)
        )

    yield (
        gr.update(value="<div style='padding: 20px; background: #fff3cd; border-radius: 8px; color: #856404;'>⚙️ Converting PDBQT to PDB...</div>", visible=True),
        gr.update(value=None, visible=False)
    )
    
    output_dir = PRANKWEB_OUTPUT_DIR
    os.makedirs(output_dir, exist_ok=True)
    
    # 2. Convert PDBQT -> PDB (P2Rank requires PDB)
    temp_pdb_filename = "prepared_for_p2rank.pdb"
    # Use ABSOLUTE paths because we are changing CWD later
    temp_pdb_path = os.path.abspath(os.path.join(output_dir, temp_pdb_filename))
    abs_output_dir = os.path.abspath(output_dir)
    
    success = convert_pdbqt_to_pdb(prepared_pdbqt_path, temp_pdb_path)
    
    if not success:
         return (
            gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Failed to convert PDBQT to PDB. Please ensure <b>OpenBabel</b> is installed.</div>", visible=True),
            gr.update(value=None, visible=False)
        )

    yield (
        gr.update(value="<div style='padding: 20px; background: #fff3cd; border-radius: 8px; color: #856404;'>⚡ Running local P2Rank prediction...</div>", visible=True),
        gr.update(value=None, visible=False)
    )

    # 3. Execute Local P2Rank
    try:
        # Command: prank predict -f "file.pdb" -o "folder"
        # We use cwd=P2RANK_HOME to simulate "cd into folder"
        cmd = [
            P2RANK_CMD,
            "predict",
            "-f", temp_pdb_path,
            "-o", abs_output_dir
        ]
        
        # Run subprocess inside the p2rank folder
        # shell=True is often required for .bat files on Windows
        process = subprocess.run(
            cmd, 
            cwd=P2RANK_HOME, 
            capture_output=True, 
            text=True, 
            shell=True 
        )

        if process.returncode != 0:
            raise Exception(f"P2Rank exited with error: {process.stderr}")

        # 4. Process Results
        # P2Rank output format is usually: [pdb_filename]_predictions.csv
        expected_csv_name = f"{os.path.splitext(temp_pdb_filename)[0]}_predictions.csv"
        csv_path = os.path.join(output_dir, expected_csv_name)

        # If exact match fails, try to find any *_predictions.csv in the output dir
        if not os.path.exists(csv_path):
             found_csvs = glob.glob(os.path.join(output_dir, "*_predictions.csv"))
             if found_csvs:
                 csv_path = found_csvs[0]
             else:
                 raise FileNotFoundError("Could not find _predictions.csv in output folder.")

        # Read and filter CSV
        df = pd.read_csv(csv_path)
        
        # Clean up column names (strip whitespace)
        df.columns = df.columns.str.strip()
        
        columns_to_drop = ['residue_ids', 'surf_atom_ids']
        df = df.drop(columns=[col for col in columns_to_drop if col in df.columns], errors='ignore')
        
        # Store CSV path globally for later use in docking
        current_pdb_info["prankweb_csv"] = csv_path
        
        yield (
            gr.update(value=f"<div style='padding: 20px; background: #d4edda; border-radius: 8px; color: #155724;'>✅ Local P2Rank prediction completed!</div>", visible=True),
            gr.update(value=df, visible=True)
        )
        
    except Exception as e:
        yield (
            gr.update(value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Error executing P2Rank: {str(e)}</div>", visible=True),
            gr.update(value=None, visible=False)
        )