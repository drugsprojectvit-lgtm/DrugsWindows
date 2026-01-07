"""
Ramachandran plot analysis module with SWISS-MODEL integration
"""

import os
import subprocess
import gradio as gr
import requests
import time
import re
from config import current_pdb_info, RAMPLOT_OUTPUT_DIR, PROTEINS_DIR

# SWISS-MODEL Configuration
SWISS_MODEL_API_TOKEN = "9e8b3ac03b851bb3834cdb311045c78021087d1d"
SWISS_MODEL_BASE_URL = "https://swissmodel.expasy.org"
SWISS_MODEL_HEADERS = {"Authorization": f"Token {SWISS_MODEL_API_TOKEN}"}


def check_remark_465(pdb_path: str) -> bool:
    """
    Check if REMARK 465 (missing residues) is present in PDB file.
    Returns True if missing residues are found, False otherwise.
    """
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('REMARK 465'):
                    # Skip header lines
                    if 'MISSING RESIDUES' in line or 'SSSEQ' in line:
                        continue
                    # If we find actual data lines, missing residues exist
                    if line.strip() and len(line.split()) > 2:
                        return True
        return False
    except Exception as e:
        print(f"Error checking REMARK 465: {e}")
        return False


def extract_favoured_percentage(csv_path: str) -> float:
    """
    Extract the Favoured percentage from CSV file.
    Looks for pattern like: Favoured: ,232,(85.294%)
    """
    try:
        with open(csv_path, 'r') as f:
            content = f.read()
            
            # Search for pattern: Favoured: ,xxx,(yy.yyy%)
            match = re.search(r'Favoured:\s*,\d+,\((\d+\.?\d*)%\)', content)
            
            if match:
                percentage = float(match.group(1))
                return percentage
            else:
                print(f"Warning: Could not find 'Favoured:' percentage in {csv_path}")
                return 0.0
                
    except FileNotFoundError:
        print(f"Warning: CSV file not found at '{csv_path}'")
        return 0.0
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return 0.0


def parse_fasta_file(fasta_path: str) -> list:
    """
    Parse FASTA file and return list of sequences.
    """
    sequences = []
    current_sequence = []
    
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                
                # If we encounter a header line
                if line.startswith('>'):
                    # If we have a sequence accumulated, save it
                    if current_sequence:
                        sequences.append("".join(current_sequence))
                        current_sequence = []
                else:
                    # Add sequence line (ignore empty lines)
                    if line:
                        current_sequence.append(line)
            
            # Don't forget the last sequence
            if current_sequence:
                sequences.append("".join(current_sequence))
        
        return sequences
    except Exception as e:
        print(f"Error parsing FASTA file: {e}")
        return []


def run_swiss_model(fasta_path: str, pdb_id: str, progress_callback=None) -> str:
    """
    Run SWISS-MODEL homology modeling on the FASTA sequence.
    Returns the path to the generated PDB file or None if failed.
    """
    if progress_callback:
        progress_callback(0.1, desc="📖 Reading FASTA file...")
    
    # Parse FASTA file
    sequences = parse_fasta_file(fasta_path)
    
    if not sequences:
        print("Error: No valid sequences found in FASTA file")
        return None
    
    # Use single sequence or multiple
    if len(sequences) == 1:
        fasta_input = sequences[0]
        print(f"Single sequence detected with {len(fasta_input)} residues")
    else:
        fasta_input = sequences
        print(f"Multiple sequences detected: {[len(s) for s in sequences]} residues each")
    
    if progress_callback:
        progress_callback(0.2, desc="🚀 Submitting job to SWISS-MODEL...")
    
    # Submit modeling job
    project_title = f"Homology_Model_{pdb_id}"
    payload = {
        "target_sequences": fasta_input,
        "project_title": project_title
    }
    
    try:
        submit_response = requests.post(
            f"{SWISS_MODEL_BASE_URL}/automodel/", 
            headers=SWISS_MODEL_HEADERS, 
            json=payload,
            timeout=60
        )
        submit_response.raise_for_status()
        
        project_id = submit_response.json().get("project_id")
        print(f"Job submitted successfully! Project ID: {project_id}")
        
    except requests.exceptions.HTTPError as e:
        print(f"Error submitting job: {e.response.status_code}")
        print(f"Details: {e.response.text}")
        return None
    except Exception as e:
        print(f"Error submitting SWISS-MODEL job: {e}")
        return None
    
    # Poll for results
    max_wait_time = 1800  # 30 minutes max
    start_time = time.time()
    poll_interval = 30  # Check every 30 seconds
    
    while True:
        elapsed = time.time() - start_time
        if elapsed > max_wait_time:
            print("SWISS-MODEL job timed out after 30 minutes")
            return None
        
        if progress_callback:
            progress_pct = 0.2 + (0.6 * (elapsed / max_wait_time))
            progress_callback(progress_pct, desc=f"⏳ Waiting for SWISS-MODEL (elapsed: {int(elapsed)}s)...")
        
        try:
            status_response = requests.get(
                f"{SWISS_MODEL_BASE_URL}/project/{project_id}/models/summary/", 
                headers=SWISS_MODEL_HEADERS,
                timeout=60
            )
            status_response.raise_for_status()
            
            status_data = status_response.json()
            job_status = status_data.get("status")
            
            print(f"Current status: {job_status}")
            
            if job_status == "COMPLETED":
                print("Modeling completed!")
                
                if progress_callback:
                    progress_callback(0.85, desc="📥 Downloading model...")
                
                # Download the PDB file
                models = status_data.get("models")
                if not models:
                    print("Job completed but no models were found")
                    return None
                
                # Get the first (and usually best) model's ID
                model_id = models[0].get("model_id")
                output_filename = os.path.join(PROTEINS_DIR, f"{pdb_id}_swiss_model.pdb")
                
                print(f"Downloading model {model_id} to {output_filename}...")
                
                pdb_response = requests.get(
                    f"{SWISS_MODEL_BASE_URL}/project/{project_id}/models/{model_id}.pdb",
                    headers=SWISS_MODEL_HEADERS,
                    timeout=60
                )
                pdb_response.raise_for_status()
                
                # Save the file
                with open(output_filename, "w") as f:
                    f.write(pdb_response.text)
                    
                print(f"Successfully saved model to: {output_filename}")
                return output_filename
                
            elif job_status == "FAILED":
                print("SWISS-MODEL job failed")
                return None
                
            elif job_status in ["RUNNING", "PENDING"]:
                print(f"Job is still running. Waiting {poll_interval} seconds...")
                time.sleep(poll_interval)
            
            else:
                print(f"Unknown status: {job_status}. Waiting...")
                time.sleep(poll_interval)
        
        except requests.exceptions.HTTPError as e:
            print(f"Error checking status: {e.response.status_code}. Retrying...")
            time.sleep(poll_interval)
        except Exception as e:
            print(f"Error during status check: {e}. Retrying...")
            time.sleep(poll_interval)


def run_ramplot(progress=gr.Progress()):
    """
    Run Ramachandran plot analysis with automatic SWISS-MODEL integration.
    If REMARK 465 is found, runs SWISS-MODEL first.
    """
    # -------------------------
    # 1️⃣ No structure case
    # -------------------------
    if not current_pdb_info["pdb_id"] or not current_pdb_info["pdb_path"]:
        return (
            gr.update(
                value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ No structure loaded. Please search for a disease/protein first.</div>",
                visible=True
            ),
            gr.update(visible=False),
            gr.update(visible=False),
            gr.update(visible=False),
            gr.update(visible=False)
        )

    pdb_id = current_pdb_info["pdb_id"]
    pdb_path = current_pdb_info["pdb_path"]
    
    # -------------------------
    # 2️⃣ Check for REMARK 465
    # -------------------------
    progress(0.05, desc="🔍 Checking for missing residues (REMARK 465)...")
    
    has_missing_residues = check_remark_465(pdb_path)
    
    if has_missing_residues:
        progress(0.1, desc="⚠️ Missing residues detected - SWISS-MODEL required...")
        
        # Check if FASTA file exists
        fasta_path = os.path.join(PROTEINS_DIR, f"{pdb_id}.fasta")
        
        if not os.path.exists(fasta_path):
            return (
                gr.update(
                    value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ FASTA file not found at {fasta_path}. Cannot proceed with SWISS-MODEL.</div>",
                    visible=True
                ),
                gr.update(visible=False),
                gr.update(visible=False),
                gr.update(visible=False),
                gr.update(visible=False)
            )
        
        # Run SWISS-MODEL
        swiss_model_path = run_swiss_model(fasta_path, pdb_id, progress)
        
        if not swiss_model_path:
            return (
                gr.update(
                    value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ SWISS-MODEL homology modeling failed. Cannot proceed.</div>",
                    visible=True
                ),
                gr.update(visible=False),
                gr.update(visible=False),
                gr.update(visible=False),
                gr.update(visible=False)
            )
        
        # Update the PDB path to use the SWISS-MODEL generated structure
        pdb_path = swiss_model_path
        current_pdb_info["pdb_path"] = swiss_model_path
        
        progress(0.9, desc="✅ SWISS-MODEL complete - proceeding with analysis...")
    else:
        progress(0.1, desc="✅ No missing residues - using original structure...")

    # -------------------------
    # 3️⃣ Run Ramachandran analysis
    # -------------------------
    progress(0.3, desc="🔬 Running Ramachandran plot analysis...")

    try:
        input_folder = PROTEINS_DIR
        output_folder = RAMPLOT_OUTPUT_DIR

        os.makedirs(input_folder, exist_ok=True)
        os.makedirs(output_folder, exist_ok=True)

        progress(0.5, desc="Executing ramplot command...")

        cmd = [
            "ramplot", "pdb",
            "-i", input_folder,
            "-o", output_folder,
            "-m", "0",
            "-r", "600",
            "-p", "png"
        ]

        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=120
        )

        progress(0.8, desc="Loading generated plots...")

        plot_dir = os.path.join(output_folder, "Plots")

        plot_files = {
            'map2d': os.path.join(plot_dir, "MapType2DAll.png"),
            'map3d': os.path.join(plot_dir, "MapType3DAll.png"),
            'std2d': os.path.join(plot_dir, "StdMapType2DGeneralGly.png"),
            'std3d': os.path.join(plot_dir, "StdMapType3DGeneral.png"),
        }

        # Safety check
        for name, path in plot_files.items():
            if not os.path.exists(path):
                raise FileNotFoundError(f"Missing plot: {path}")

        progress(1.0, desc="✅ Complete!")

        # -------------------------
        # 4️⃣ Final successful UI update
        # -------------------------
        success_msg = "✅ Ramachandran plot analysis completed!"
        if has_missing_residues:
            success_msg += "<br>🔧 Used SWISS-MODEL homology model (original had missing residues)"
        
        return (
            gr.update(
                value=f"<div style='padding: 20px; background: #d4edda; border-radius: 8px; color: #155724;'>{success_msg}</div>",
                visible=True
            ),
            gr.update(value=plot_files['map2d'], visible=True),
            gr.update(value=plot_files['map3d'], visible=True),
            gr.update(value=plot_files['std2d'], visible=True),
            gr.update(value=plot_files['std3d'], visible=True)
        )

    except subprocess.CalledProcessError as e:
        return (
            gr.update(
                value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>⚠️ Ramachandran analysis failed:<br><pre>{e.stderr.decode()}</pre></div>",
                visible=True
            ),
            gr.update(visible=False),
            gr.update(visible=False),
            gr.update(visible=False),
            gr.update(visible=False)
        )

    except Exception as e:
        return (
            gr.update(
                value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Error: {str(e)}</div>",
                visible=True
            ),
            gr.update(visible=False),
            gr.update(visible=False),
            gr.update(visible=False),
            gr.update(visible=False)
        )