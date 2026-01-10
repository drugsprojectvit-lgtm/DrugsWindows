# batch.py
"""
Batch Pipeline for Protein Structure Analysis
Processes multiple proteins through the complete workflow and saves all outputs.
UPDATED:
- Compatible with Local P2Rank (prankweb.py).
- Compatible with Clean Docking workflow (docking.py).
- Compatible with robust ADMET analysis (admet_analysis.py).
- ADDED: FASTA download in Step 1 to support SWISS-MODEL fallback in Step 2.
"""

import os
import json
import shutil
from pathlib import Path
from datetime import datetime
import pandas as pd
import traceback
import time
import glob
import base64
import re
import requests # Added for FASTA download

# Import modules from your app
from config import current_pdb_info, PROTEINS_DIR, DOCKING_RESULTS_DIR, RAMPLOT_OUTPUT_DIR
from ramachandran import run_ramplot
from prankweb import run_prankweb_prediction
from protein_prep import prepare_protein_meeko
from docking import run_molecular_docking
from admet_analysis import run_admet_prediction
from utils import map_disease_to_protein, find_best_pdb_structure


# --- VISUALIZATION FUNCTION ---
def show_structure(protein_text: str, ligand_text: str = None, pdb_id: str = "Docking Result", protein_name: str = "") -> str:
    """Robust 3D visualization supporting Overlay."""
    prot_json = json.dumps(protein_text)
    lig_json = json.dumps(ligand_text) if ligand_text else "null"
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <style>
            body {{ margin: 0; padding: 0; overflow: hidden; background-color: white; }}
            #container {{ width: 100vw; height: 100vh; position: relative; }}
            #error-log {{ 
                display: none; position: absolute; top: 10px; left: 10px; 
                background: rgba(255,0,0,0.8); color: white; padding: 10px; 
                z-index: 999; font-family: sans-serif; border-radius: 5px;
            }}
            .legend {{
                position: absolute; bottom: 10px; right: 10px;
                background: rgba(255, 255, 255, 0.9); padding: 8px 12px;
                border-radius: 6px; font-family: sans-serif; font-size: 12px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.2); z-index: 100;
                pointer-events: none;
            }}
            .color-box {{ 
                display: inline-block; width: 10px; height: 10px; 
                margin-right: 6px; border-radius: 50%; 
            }}
            .info-overlay {{
                position: absolute; top: 10px; left: 10px;
                background: rgba(255, 255, 255, 0.9); padding: 8px 12px;
                border-radius: 6px; font-family: sans-serif; font-size: 14px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.2); z-index: 90;
            }}
        </style>
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    </head>
    <body>
        <div id="error-log"></div>
        <div class="info-overlay">
            <strong>{pdb_id}</strong><br>
            <span style="font-size:12px; color:#555">{protein_name}</span>
        </div>
        <div class="legend">
            <div><span class="color-box" style="background: linear-gradient(90deg, blue, green, red);"></span>Protein</div>
            <div style="margin-top:4px"><span class="color-box" style="background: magenta;"></span>Ligand</div>
        </div>
        <div id="container"></div>
        <script>
            function logError(msg) {{
                var el = document.getElementById('error-log');
                el.style.display = 'block';
                el.innerHTML += "Error: " + msg + "<br>";
                console.error(msg);
            }}
            window.onload = function() {{
                try {{
                    if (typeof $3Dmol === 'undefined') throw new Error("3Dmol.js failed to load.");
                    var element = document.getElementById('container');
                    var config = {{ backgroundColor: 'white' }};
                    var viewer = $3Dmol.createViewer(element, config);

                    var proteinData = {prot_json};
                    if (!proteinData) throw new Error("Protein data is empty.");
                    
                    viewer.addModel(proteinData, "pdb");
                    viewer.setStyle({{model: 0}}, {{cartoon: {{color: 'spectrum'}}}});
                    viewer.addStyle({{model: 0, hetflag: true}}, {{stick: {{radius: 0.1, color: 'lightgray'}}}});

                    var ligandData = {lig_json};
                    if (ligandData) {{
                        viewer.addModel(ligandData, "pdb");
                        viewer.setStyle({{model: 1}}, {{stick: {{colorscheme: 'magentaCarbon', radius: 0.4}}}});
                    }} else {{
                        viewer.addStyle({{resn: ["UNL", "LIG", "DRG", "UNK"]}}, {{stick: {{colorscheme: 'magentaCarbon', radius: 0.4}}}});
                    }}

                    viewer.zoomTo();
                    viewer.render();
                }} catch(e) {{ logError(e.message); }}
            }};
        </script>
    </body>
    </html>
    """
    b64 = base64.b64encode(html_content.encode()).decode()
    iframe = f'<iframe src="data:text/html;base64,{b64}" width="100%" height="500" frameborder="0" style="border: 1px solid #ccc; border-radius: 8px;"></iframe>'
    return iframe


class ProteinPipelineBatch:
    """Handles batch processing of proteins through the complete pipeline."""
    
    def __init__(self, output_base_dir="batch_results"):
        """Initialize batch processor with output directory."""
        self.output_base_dir = Path(output_base_dir)
        self.output_base_dir.mkdir(exist_ok=True)
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
    def create_protein_folders(self, protein_name):
        """Create folder structure for a protein."""
        safe_name = "".join(c if c.isalnum() or c in (' ', '-', '_') else '_' 
                           for c in protein_name).strip()
        protein_dir = self.output_base_dir / safe_name
        protein_dir.mkdir(exist_ok=True)
        
        steps = ["01_structure_search", "02_ramachandran_analysis", "03_protein_preparation",
                 "04_binding_site_prediction", "05_molecular_docking", "06_admet_analysis"]
        step_dirs = {}
        for step in steps:
            step_dir = protein_dir / step
            step_dir.mkdir(exist_ok=True)
            step_dirs[step] = step_dir
        return protein_dir, step_dirs
    
    def save_3d_viewer_as_image(self, html_content, output_path, title=""):
        """Convert 3D viewer HTML to 2D image using screenshot method with WebGL support."""
        if not html_content: return False

        if isinstance(html_content, str) and "data:text/html;base64," in html_content:
            try:
                # Decode the base64 iframe content to get the actual HTML
                match = re.search(r'src="data:text/html;base64,([^"]+)"', html_content)
                if match:
                    b64_str = match.group(1)
                    decoded_bytes = base64.b64decode(b64_str)
                    html_content = decoded_bytes.decode('utf-8')
            except Exception as e:
                print(f"⚠️ Failed to decode base64 iframe: {e}")
        
        try:
            from selenium import webdriver
            from selenium.webdriver.chrome.options import Options
            from selenium.webdriver.common.by import By
            
            chrome_options = Options()
            chrome_options.add_argument("--headless=new") 
            chrome_options.add_argument("--use-gl=swiftshader") 
            chrome_options.add_argument("--enable-webgl")
            chrome_options.add_argument("--ignore-gpu-blocklist")
            chrome_options.add_argument("--no-sandbox")
            chrome_options.add_argument("--disable-dev-shm-usage")
            chrome_options.add_argument("--window-size=1920,1080")
            
            driver = webdriver.Chrome(options=chrome_options)
            
            style_injection = "<style>body, html, .viewer-container, #viewport, #container { width: 100%; height: 100%; margin: 0; padding: 0; overflow: hidden; }</style>"
            if "</head>" in html_content:
                html_content = html_content.replace("</head>", f"{style_injection}</head>")
            else:
                html_content = style_injection + html_content

            temp_html = output_path.parent / "temp_viewer.html"
            with open(temp_html, 'w', encoding='utf-8') as f: f.write(html_content)
            
            driver.get(f"file://{temp_html.absolute()}")
            time.sleep(8)  # Wait for 3Dmol to render
            
            try: driver.execute_script("if(typeof viewer !== 'undefined') { viewer.render(); }")
            except: pass

            # Try to screenshot the canvas specifically, else full page
            try:
                canvas = driver.find_element(By.TAG_NAME, "canvas")
                canvas.screenshot(str(output_path))
            except:
                driver.save_screenshot(str(output_path))
                
            driver.quit()
            if temp_html.exists(): temp_html.unlink()
            return True
        except Exception as e:
            print(f"⚠️ Screenshot failed: {e}")
            # Save HTML as fallback if screenshot fails
            html_path = output_path.with_suffix('.html')
            with open(html_path, 'w', encoding='utf-8') as f: f.write(html_content)
            return False
    
    # ... [Steps 1-3] ...

    def process_structure_search(self, protein_input, step_dir):
        print(f"\n{'='*60}")
        print(f"STEP 1: Structure Search - {protein_input}")
        print(f"{'='*60}")
        results = {"status": "failed", "protein_input": protein_input}
        try:
            protein_name = map_disease_to_protein(protein_input) or protein_input.strip()
            results["protein_name"] = protein_name
            result = find_best_pdb_structure(protein_name, max_check=100)
            if not result:
                results["error"] = "No suitable PDB structure found"
                print(f"❌ No structure found")
                return results
            
            pdb_id, pdb_path = result
            # Ensure Absolute Path
            pdb_path = os.path.abspath(pdb_path)
            results["pdb_id"], results["pdb_path"] = pdb_id, pdb_path
            
            current_pdb_info.update({"pdb_id": pdb_id, "pdb_path": pdb_path, "prepared_pdbqt": None, "docking_results": None, "prankweb_csv": None})
            
            dest_pdb = step_dir / f"{pdb_id}.pdb"
            shutil.copy2(pdb_path, dest_pdb)
            
            # --- NEW: Download FASTA for SWISS-MODEL Support ---
            try:
                print(f"  ⬇️  Fetching FASTA for {pdb_id} (Required for potential Homology Modeling)...")
                fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
                response = requests.get(fasta_url, timeout=10)
                
                if response.status_code == 200:
                    fasta_content = response.text
                    
                    # 1. Save to PROTEINS_DIR (Required by ramachandran.py logic)
                    os.makedirs(PROTEINS_DIR, exist_ok=True)
                    main_fasta_path = os.path.join(PROTEINS_DIR, f"{pdb_id}.fasta")
                    with open(main_fasta_path, "w") as f:
                        f.write(fasta_content)
                        
                    # 2. Save to Batch Step Directory (For record keeping)
                    batch_fasta_path = step_dir / f"{pdb_id}.fasta"
                    with open(batch_fasta_path, "w") as f:
                        f.write(fasta_content)
                        
                    print(f"  ✅ FASTA saved to {main_fasta_path}")
                else:
                    print(f"  ⚠️ Failed to download FASTA (Status: {response.status_code})")
            except Exception as fasta_err:
                print(f"  ⚠️ Error downloading FASTA: {fasta_err}")
            # ----------------------------------------------------

            with open(pdb_path, 'r') as f: pdb_content = f.read()
            structure_html = show_structure(protein_text=pdb_content, pdb_id=pdb_id, protein_name=protein_name)
            self.save_3d_viewer_as_image(structure_html, step_dir / f"{pdb_id}_structure.png", f"{protein_name}")
            results["status"] = "success"
            print(f"✅ Structure loaded: {pdb_id}")
        except Exception as e:
            results["error"] = str(e)
            print(f"❌ Error: {e}")
        return results

    def process_ramachandran(self, step_dir):
        print(f"\nSTEP 2: Ramachandran Analysis")
        run_ramplot() 
        
        csv_files_copied = 0
        img_files_copied = 0
        
        if os.path.exists(RAMPLOT_OUTPUT_DIR):
            for csv_file in glob.glob(os.path.join(RAMPLOT_OUTPUT_DIR, "*.csv")):
                try:
                    shutil.copy2(csv_file, step_dir)
                    csv_files_copied += 1
                except Exception as e:
                    print(f"Warning: Could not copy csv {csv_file}: {e}")
            
            plots_subdir = os.path.join(RAMPLOT_OUTPUT_DIR, "Plots")
            target_images = ["MapType2DAll.png", "MapType3DAll.png", "StdMapType2DGeneralGly.png", "StdMapType3DGeneral.png"]
            
            if os.path.exists(plots_subdir):
                for img_name in target_images:
                    src_img = os.path.join(plots_subdir, img_name)
                    if os.path.exists(src_img):
                        try:
                            shutil.copy2(src_img, step_dir)
                            img_files_copied += 1
                        except Exception as e:
                            print(f"Warning: Could not copy image {img_name}: {e}")
            else:
                print(f"Warning: 'Plots' subdirectory not found in {RAMPLOT_OUTPUT_DIR}")
        
        print(f"✅ Ramachandran analysis complete. Copied {csv_files_copied} CSVs and {img_files_copied} images.")
        return {"status": "success", "csv_files_count": csv_files_copied, "img_files_count": img_files_copied}

    def process_protein_preparation(self, step_dir):
        print(f"\nSTEP 3: Protein Preparation")
        generator = prepare_protein_meeko()
        final_output = None
        for output in generator: final_output = output
        
        results = {"status": "failed"}
        if final_output:
            pdbqt_update = final_output[2] if len(final_output) > 2 else None
            # Handle potential Gradio update object or direct dict
            if hasattr(pdbqt_update, 'value'):
                pdbqt_path = pdbqt_update.value
            elif isinstance(pdbqt_update, dict):
                pdbqt_path = pdbqt_update.get('value')
            else:
                pdbqt_path = pdbqt_update

            if pdbqt_path and os.path.exists(pdbqt_path):
                dest_pdbqt = step_dir / Path(pdbqt_path).name
                shutil.copy2(pdbqt_path, dest_pdbqt)
                current_pdb_info["prepared_pdbqt"] = str(dest_pdbqt)
                
                try:
                    with open(pdbqt_path, 'r') as f: prepared_text = f.read()
                    html_view = show_structure(protein_text=prepared_text, pdb_id="Prepared Structure", protein_name="Auto-Prepared")
                    img_path = step_dir / "prepared_protein_structure.png"
                    self.save_3d_viewer_as_image(html_view, img_path, "Prepared Structure")
                    print(f"  ✓ Saved prepared structure image")
                except Exception as img_err:
                    print(f"  ⚠️ Could not save prepared structure image: {img_err}")

                results["status"] = "success"
                print(f"✅ Protein prepared")
        return results

    def process_binding_sites(self, step_dir):
        """Step 4: Binding Site Prediction (Local P2Rank)"""
        print(f"\nSTEP 4: Binding Site Prediction")
        if not current_pdb_info.get("prepared_pdbqt"): return {"status": "skipped"}
        
        try:
            print(f"  Running Local PrankWeb (P2Rank)...")
            generator = run_prankweb_prediction()
            final_output = None
            
            # Iterate through generator to trigger execution
            for output in generator: final_output = output
            
            results = {"status": "failed"}
            
            if final_output:
                # In your new prankweb.py, output is (status_update, df_update)
                df_update = final_output[1] if len(final_output) > 1 else None
                
                # Extract DataFrame from Gradio update object
                df = None
                if hasattr(df_update, 'value'):
                    df = df_update.value
                elif isinstance(df_update, dict):
                    df = df_update.get('value')
                else:
                    df = df_update

                if df is not None and not df.empty:
                    csv_path = step_dir / "binding_sites.csv"
                    df.to_csv(csv_path, index=False)
                    current_pdb_info["prankweb_csv"] = str(csv_path)
                    results["status"] = "success"
                    print(f"✅ {len(df)} pockets found")
                    return results # SUCCESS
                else:
                     print("  ⚠️ Output DataFrame is empty.")
            else:
                  print("  ⚠️ No final output received from P2Rank.")

        except Exception as e:
            print(f"  ❌ Error: {e}")
            traceback.print_exc()

        return {"status": "failed", "error": "P2Rank execution failed"}

    def process_docking(self, step_dir):
        """Step 5: Molecular docking (ALL RESULTS - NO FILTER)"""
        print(f"\n{'='*60}")
        print(f"STEP 5: Molecular Docking")
        print(f"{'='*60}")
        
        results = {"status": "failed"}
        
        try:
            if not current_pdb_info.get("prepared_pdbqt"): return results
            if not current_pdb_info.get("prankweb_csv"): return results

            generator = run_molecular_docking()
            final_output = None
            print("  ⚙️  Running Vina Docking...")
            for output in generator: final_output = output
            
            if not final_output: return results
            
            # Extract summary dataframe
            df_update = final_output[1] if len(final_output) > 1 else None
            summary_df = None
            if hasattr(df_update, 'value'):
                summary_df = df_update.value
            elif isinstance(df_update, dict):
                summary_df = df_update.get('value')
            else:
                summary_df = df_update
            
            if summary_df is not None and not summary_df.empty:
                csv_path = step_dir / "docking_summary.csv"
                summary_df.to_csv(csv_path, index=False)
                results["csv_file"] = str(csv_path)
                print(f"  ✓ Saved: docking_summary.csv")
                
                summary_df['binding_energy'] = pd.to_numeric(summary_df['binding_energy'], errors='coerce')
                
                # Process ALL results
                results_to_process = summary_df 
                print(f"  ℹ️  Processing all {len(results_to_process)} docked poses...")

                # Pre-load Protein Structure for Visualization
                protein_path = current_pdb_info.get("pdb_path")
                protein_text = ""
                if protein_path and os.path.exists(protein_path):
                    with open(protein_path, 'r') as pf: protein_text = pf.read()
                else:
                    print("  ⚠️ Protein structure not found for visualization.")

                for idx, row in results_to_process.iterrows():
                    ligand_name = str(row['ligand']).strip()
                    pocket_name = str(row['pocket']).strip()
                    energy = float(row['binding_energy'])
                    pose_num = int(row['pose_number'])
                    pdb_file_path = row['pdb_file'] # Ligand Only PDB
                    
                    safe_ligand = "".join(c if c.isalnum() else '_' for c in ligand_name)
                    energy_str = f"{energy:.1f}kcal"
                    base_filename = f"{safe_ligand}_{energy_str}_pose_{pose_num}"

                    # 1. Save 2D Interaction Image (Generated by pandamap in docking.py)
                    if 'interaction_image' in row and row['interaction_image'] != "N/A":
                        img_src = row['interaction_image']
                        if os.path.exists(img_src):
                            img_dest = step_dir / f"{base_filename}_2D_interaction.png"
                            shutil.copy2(img_src, img_dest)
                            print(f"  ✓ Saved 2D Image: {img_dest.name}")

                    # 2. SAVE COMPLEX PDB
                    output_dir_pdb = os.path.join(DOCKING_RESULTS_DIR, "pdb")
                    complex_src = os.path.join(output_dir_pdb, f"{ligand_name}_{pocket_name}_complex.pdb")
                    if os.path.exists(complex_src):
                        complex_dest = step_dir / f"{base_filename}_complex.pdb"
                        shutil.copy2(complex_src, complex_dest)
                        print(f"  ✓ Saved Complex PDB: {complex_dest.name}")
                    else:
                        print(f"  ⚠️ Complex PDB not found: {complex_src}")

                    # 3. Generate 3D Screenshot DIRECTLY
                    # We extract the specific model text manually to avoid import errors
                    ligand_text_pose = ""
                    if os.path.exists(pdb_file_path):
                        with open(pdb_file_path, 'r') as lf:
                            lines = lf.readlines()
                            model_lines = []
                            in_model = False
                            current_model = -1
                            found_pose = False
                            for line in lines:
                                if line.startswith("MODEL"):
                                    try: current_model = int(line.split()[1])
                                    except: pass
                                    if current_model == pose_num:
                                        in_model = True
                                        found_pose = True
                                if in_model: model_lines.append(line)
                                if line.startswith("ENDMDL") and in_model: break
                            
                            ligand_text_pose = "".join(model_lines) if found_pose else "".join(lines)
                    
                    if protein_text and ligand_text_pose:
                        # Generate HTML locally
                        viewer_html = show_structure(
                            protein_text=protein_text, 
                            ligand_text=ligand_text_pose, 
                            pdb_id="Docking Result", 
                            protein_name=f"Pose {pose_num} ({energy} kcal/mol)"
                        )
                        
                        screenshot_path = step_dir / f"{base_filename}_3D_view.png"
                        success = self.save_3d_viewer_as_image(viewer_html, screenshot_path, f"{ligand_name} {energy}")
                        if success:
                            print(f"  ✓ Saved 3D Screenshot: {screenshot_path.name}")
                        else:
                             print(f"  ⚠️ Failed 3D Screenshot: {base_filename}")
                    else:
                        print(f"  ⚠️ Missing protein or ligand text for 3D view: {base_filename}")

                results["status"] = "success"
                results["num_poses"] = len(results_to_process)
                print(f"✅ Docking processing complete.")
            else:
                results["error"] = "Docking failed or no results"
                print(f"❌ Docking failed - no poses found")
            
        except Exception as e:
            results["error"] = str(e)
            results["traceback"] = traceback.format_exc()
            print(f"❌ Error in docking: {e}")
            traceback.print_exc()
        
        return results

    def process_admet(self, step_dir):
        """Step 6: ADMET Analysis"""
        print(f"\nSTEP 6: ADMET Analysis")
        # run_admet_prediction now returns a tuple: (msg, df, csv_path)
        result = run_admet_prediction()
        
        if result:
            msg, df, csv_path = result
            if csv_path and os.path.exists(csv_path):
                shutil.copy2(csv_path, step_dir / "admet_results.csv")
                print(f"✅ ADMET results saved.")
                return {"status": "success", "message": msg}
            else:
                print(f"⚠️ ADMET analysis ran, but CSV path not found.")
        
        return {"status": "failed"}
    
    def process_single_protein(self, protein_input):
        print(f"\n{'#'*80}")
        print(f"# Processing: {protein_input}")
        print(f"{'#'*80}")
        
        protein_dir, step_dirs = self.create_protein_folders(protein_input)
        pipeline_results = {"protein_input": protein_input, "start_time": datetime.now().isoformat(), "steps": {}}
        
        pipeline_results["steps"]["structure_search"] = self.process_structure_search(protein_input, step_dirs["01_structure_search"])
        if pipeline_results["steps"]["structure_search"]["status"] != "success":
            pipeline_results["pipeline_status"] = "failed_at_step_1"
            self.save_pipeline_summary(protein_dir, pipeline_results)
            return pipeline_results
        
        pipeline_results["steps"]["ramachandran"] = self.process_ramachandran(step_dirs["02_ramachandran_analysis"])
        pipeline_results["steps"]["protein_preparation"] = self.process_protein_preparation(step_dirs["03_protein_preparation"])
        if pipeline_results["steps"]["protein_preparation"]["status"] != "success":
            pipeline_results["pipeline_status"] = "failed_at_step_3"
            self.save_pipeline_summary(protein_dir, pipeline_results)
            return pipeline_results
        
        pipeline_results["steps"]["binding_sites"] = self.process_binding_sites(step_dirs["04_binding_site_prediction"])
        pipeline_results["steps"]["docking"] = self.process_docking(step_dirs["05_molecular_docking"])
        pipeline_results["steps"]["admet"] = self.process_admet(step_dirs["06_admet_analysis"])
        
        pipeline_results["pipeline_status"] = "completed"
        pipeline_results["end_time"] = datetime.now().isoformat()
        
        self.save_pipeline_summary(protein_dir, pipeline_results)
        print(f"\n{'='*60}\n✅ Pipeline complete for {protein_input}\nResults saved to: {protein_dir}\n{'='*60}\n")
        return pipeline_results
    
    def save_pipeline_summary(self, protein_dir, results):
        with open(protein_dir / "pipeline_summary.json", 'w') as f: json.dump(results, f, indent=2)
    
    def run_batch(self, protein_list):
        print(f"\n{'#'*80}\n# BATCH PROCESSING: {len(protein_list)} proteins\n{'#'*80}\n")
        batch_results = {"proteins": {}}
        for i, protein in enumerate(protein_list, 1):
            print(f"[{i}/{len(protein_list)}]")
            try: batch_results["proteins"][protein] = self.process_single_protein(protein)
            except Exception as e:
                print(f"❌ Fatal error: {e}")
                batch_results["proteins"][protein] = {"status": "fatal_error", "error": str(e)}
        return batch_results
    
    def generate_batch_csv_summary(self, batch_results):
        summary_data = []
        for p, r in batch_results["proteins"].items():
            summary_data.append({"protein": p, "status": r.get("pipeline_status", "unknown")})
        pd.DataFrame(summary_data).to_csv(self.output_base_dir / f"batch_summary_{self.timestamp}.csv", index=False)


if __name__ == "__main__":
    protein_list = [
    "AKT1",
    "KRAS"
]

    batch_processor = ProteinPipelineBatch(output_base_dir="FDA_drugs")
    results = batch_processor.run_batch(protein_list)