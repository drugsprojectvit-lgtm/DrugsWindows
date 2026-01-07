"""
Batch Pipeline for Protein Structure Analysis
Processes multiple proteins through the complete workflow and saves all outputs.
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

# Import modules from your app
from config import current_pdb_info, PROTEINS_DIR, DOCKING_RESULTS_DIR
from ramachandran import run_ramplot
from prankweb import run_prankweb_prediction
from protein_prep import prepare_protein_meeko
from docking import run_molecular_docking, display_docked_structure
from admet_analysis import run_admet_prediction
from utils import map_disease_to_protein, find_best_pdb_structure
from visualization import show_structure


class ProteinPipelineBatch:
    """Handles batch processing of proteins through the complete pipeline."""
    
    def __init__(self, output_base_dir="batch_results"):
        """Initialize batch processor with output directory."""
        self.output_base_dir = Path(output_base_dir)
        self.output_base_dir.mkdir(exist_ok=True)
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
    def create_protein_folders(self, protein_name):
        """Create folder structure for a protein."""
        # Sanitize protein name for folder name
        safe_name = "".join(c if c.isalnum() or c in (' ', '-', '_') else '_' 
                           for c in protein_name).strip()
        
        protein_dir = self.output_base_dir / safe_name
        protein_dir.mkdir(exist_ok=True)
        
        # Create subfolders for each step
        steps = [
            "01_structure_search",
            "02_ramachandran_analysis",
            "03_protein_preparation",
            "04_binding_site_prediction",
            "05_molecular_docking",
            "06_admet_analysis"
        ]
        
        step_dirs = {}
        for step in steps:
            step_dir = protein_dir / step
            step_dir.mkdir(exist_ok=True)
            step_dirs[step] = step_dir
        
        return protein_dir, step_dirs
    
    def save_3d_viewer_as_image(self, html_content, output_path, title=""):
        """Convert 3D viewer HTML to 2D image using screenshot method with WebGL support."""
        try:
            from selenium import webdriver
            from selenium.webdriver.chrome.options import Options
            from selenium.webdriver.common.by import By
            from selenium.webdriver.support.ui import WebDriverWait
            from selenium.webdriver.support import expected_conditions as EC
            
            chrome_options = Options()
            
            # --- CRITICAL FIXES FOR 3D RENDERING ---
            # Use 'new' headless mode (better rendering)
            chrome_options.add_argument("--headless=new") 
            
            # Force software rendering for WebGL (Fixes blank canvas on servers/headless)
            chrome_options.add_argument("--use-gl=swiftshader") 
            chrome_options.add_argument("--enable-webgl")
            chrome_options.add_argument("--ignore-gpu-blocklist")
            
            # Standard stability flags
            chrome_options.add_argument("--no-sandbox")
            chrome_options.add_argument("--disable-dev-shm-usage")
            chrome_options.add_argument("--hide-scrollbars")
            chrome_options.add_argument("--force-device-scale-factor=1")
            
            # Set a fixed window size
            chrome_options.add_argument("--window-size=1920,1080")
            
            driver = webdriver.Chrome(options=chrome_options)
            
            # Inject CSS to ensure the viewer fills the viewport
            # (Often 3D viewers collapse to height:0 without explicit CSS)
            style_injection = """
            <style>
                body, html, .viewer-container, #viewport { 
                    width: 100%; height: 100%; margin: 0; padding: 0; overflow: hidden; 
                }
            </style>
            """
            # Insert style before closing head or at start of body
            if "</head>" in html_content:
                html_content = html_content.replace("</head>", f"{style_injection}</head>")
            else:
                html_content = style_injection + html_content

            # Create temporary HTML file
            temp_html = output_path.parent / "temp_viewer.html"
            with open(temp_html, 'w', encoding='utf-8') as f:
                f.write(html_content)
            
            # Load file
            driver.get(f"file://{temp_html.absolute()}")
            
            # Wait slightly longer for WebGL to initialize
            time.sleep(8) 
            
            # Optional: Execute JS to force a re-render if using Py3Dmol
            try:
                driver.execute_script("if(typeof viewer !== 'undefined') { viewer.render(); }")
            except:
                pass

            # Find the canvas element (usually where the 3D structure is)
            # If standard screenshot is blank, we try to grab just the canvas
            try:
                canvas = driver.find_element(By.TAG_NAME, "canvas")
                canvas.screenshot(str(output_path))
            except:
                # Fallback to full page screenshot
                driver.save_screenshot(str(output_path))
                
            driver.quit()
            
            # Remove temp file
            if temp_html.exists():
                temp_html.unlink()
            
            return True
            
        except Exception as e:
            print(f"⚠️ Screenshot failed: {e}")
            # Fallback: save HTML instead
            html_path = output_path.with_suffix('.html')
            with open(html_path, 'w', encoding='utf-8') as f:
                f.write(html_content)
            return False
    
    def process_structure_search(self, protein_input, step_dir):
        """Step 1: Search for protein structure."""
        print(f"\n{'='*60}")
        print(f"STEP 1: Structure Search - {protein_input}")
        print(f"{'='*60}")
        
        results = {
            "status": "failed",
            "protein_input": protein_input,
            "protein_name": None,
            "pdb_id": None,
            "pdb_path": None
        }
        
        try:
            # Map disease to protein if needed
            protein_name = map_disease_to_protein(protein_input)
            if not protein_name:
                protein_name = protein_input.strip()
            
            results["protein_name"] = protein_name
            
            # Find best PDB structure
            result = find_best_pdb_structure(protein_name, max_check=100)
            
            if not result:
                results["error"] = "No suitable PDB structure found"
                print(f"❌ No structure found for {protein_name}")
                return results
            
            pdb_id, pdb_path = result
            results["pdb_id"] = pdb_id
            results["pdb_path"] = pdb_path
            
            # Update global config 
            current_pdb_info.update({
                "pdb_id": pdb_id,
                "pdb_path": pdb_path,
                "prepared_pdbqt": None,
                "docking_results": None,
                "prankweb_csv": None
            })
            
            # Copy PDB file to output
            dest_pdb = step_dir / f"{pdb_id}.pdb"
            shutil.copy2(pdb_path, dest_pdb)
            
            # Read and save structure content
            with open(pdb_path, 'r') as f:
                pdb_content = f.read()
            
            # Generate 3D viewer HTML
            structure_html = show_structure(pdb_content, pdb_id, protein_name)
            
            # Save viewer as image (or HTML if screenshot fails)
            viewer_img = step_dir / f"{pdb_id}_structure.png"
            self.save_3d_viewer_as_image(structure_html, viewer_img, f"{protein_name} Structure")
            
            # Save metadata
            metadata = {
                "protein_input": protein_input,
                "protein_name": protein_name,
                "pdb_id": pdb_id,
                "pdb_file": str(dest_pdb),
                "timestamp": datetime.now().isoformat()
            }
            
            with open(step_dir / "metadata.json", 'w') as f:
                json.dump(metadata, f, indent=2)
            
            results["status"] = "success"
            print(f"✅ Structure loaded: {pdb_id}")
            
        except Exception as e:
            results["error"] = str(e)
            results["traceback"] = traceback.format_exc()
            print(f"❌ Error in structure search: {e}")
        
        return results
    
    def process_ramachandran(self, step_dir):
        """Step 2: Ramachandran analysis (With Fallback)."""
        print(f"\n{'='*60}")
        print(f"STEP 2: Ramachandran Analysis")
        print(f"{'='*60}")
        
        results = {"status": "failed"}
        
        # 1. CAPTURE ORIGINAL PDB (Safety Net)
        original_pdb = current_pdb_info.get("pdb_path")
        
        try:
            # Run ramachandran analysis
            outputs = run_ramplot()
            
            # Extract values from Gradio update objects
            status_update = outputs[0]
            plot_updates = outputs[1:5]
            
            status_html = status_update.get('value', '') if isinstance(status_update, dict) else str(status_update)
            
            # Save plots 
            plot_files = []
            plot_names = ["ramachandran_map_2d.png", "ramachandran_map_3d.png", 
                          "ramachandran_std_2d.png", "ramachandran_std_3d.png"]
            
            for i, (plot_update, plot_name) in enumerate(zip(plot_updates, plot_names)):
                # Extract the actual file path from Gradio update
                if isinstance(plot_update, dict):
                    plot_path = plot_update.get('value')
                else:
                    plot_path = plot_update
                
                if plot_path and os.path.exists(plot_path):
                    dest_path = step_dir / plot_name
                    shutil.copy2(plot_path, dest_path)
                    plot_files.append(str(dest_path))
            
            if plot_files:
                results["status"] = "success"
                results["plot_files"] = plot_files
                results["status_message"] = status_html
                print(f"✅ Ramachandran analysis complete: {len(plot_files)} plots saved")
            else:
                # 2. FALLBACK LOGIC TRIGGERED
                results["error"] = "No plots were generated"
                print(f"⚠️  Homology modeling finished but returned no plots.")
                
                if original_pdb and os.path.exists(original_pdb):
                    # Revert config to use original PDB
                    current_pdb_info["pdb_path"] = original_pdb
                    
                    results["status"] = "fallback"
                    results["message"] = "Used original PDB fallback"
                    print(f"ℹ️  FALLBACK: Reverting to original Crystal Structure.")
                    print(f"   Path: {original_pdb}")
                    print(f"   The pipeline will continue using the original PDB file.")
                else:
                    print("❌ Critical Error: Original PDB file is missing.")
            
        except Exception as e:
            # 3. EXCEPTION FALLBACK
            print(f"⚠️ Error during Ramachandran step: {e}")
            if original_pdb and os.path.exists(original_pdb):
                current_pdb_info["pdb_path"] = original_pdb
                results["status"] = "fallback_error"
                print(f"ℹ️  FALLBACK (Exception): Reverting to original PDB structure.")
            else:
                results["error"] = str(e)
                results["traceback"] = traceback.format_exc()
        
        return results
    
    def process_protein_preparation(self, step_dir):
        """Step 3: Protein preparation."""
        print(f"\n{'='*60}")
        print(f"STEP 3: Protein Preparation")
        print(f"{'='*60}")
        
        results = {"status": "failed"}
        
        try:
            # Run protein preparation (CONSUME GENERATOR)
            generator = prepare_protein_meeko()
            final_output = None
            
            print("  ⚙️  Running Meeko preparation...", end="", flush=True)
            for output in generator:
                final_output = output
                print(".", end="", flush=True)
            print(" Done.")
            
            if not final_output:
                results["error"] = "No output from protein preparation"
                print("\n❌ Protein preparation returned no output")
                return results

            # Extract values from final Gradio update objects
            status_update = final_output[0]
            viewer_update = final_output[1] if len(final_output) > 1 else None
            pdbqt_update = final_output[2] if len(final_output) > 2 else None
            
            status_html = status_update.get('value', '') if isinstance(status_update, dict) else str(status_update)
            viewer_html = viewer_update.get('value') if isinstance(viewer_update, dict) else viewer_update
            
            # Extract PDBQT path
            if isinstance(pdbqt_update, dict):
                pdbqt_path = pdbqt_update.get('value')
            else:
                pdbqt_path = pdbqt_update
            
            if pdbqt_path and os.path.exists(pdbqt_path):
                # Copy PDBQT file to batch directory
                dest_pdbqt = step_dir / Path(pdbqt_path).name
                shutil.copy2(pdbqt_path, dest_pdbqt)
                results["pdbqt_file"] = str(dest_pdbqt)
                
                # CRITICAL: Update current_pdb_info to point to THIS batch file
                current_pdb_info["prepared_pdbqt"] = str(dest_pdbqt)
                
                print(f"  ✓ Saved PDBQT: {dest_pdbqt.name}")
                
                # Save viewer
                if viewer_html:
                    viewer_img = step_dir / "prepared_structure.png"
                    self.save_3d_viewer_as_image(viewer_html, viewer_img, "Prepared Structure")
                
                results["status"] = "success"
                results["status_message"] = status_html
                print(f"✅ Protein preparation complete")
            else:
                results["error"] = "Failed to generate PDBQT file"
                results["status_message"] = status_html
                print(f"❌ Protein preparation failed - no PDBQT file generated")
            
        except Exception as e:
            results["error"] = str(e)
            results["traceback"] = traceback.format_exc()
            print(f"❌ Error in protein preparation: {e}")
        
        return results
    
    def process_binding_sites(self, step_dir):
        """Step 4: Binding site prediction."""
        print(f"\n{'='*60}")
        print(f"STEP 4: Binding Site Prediction")
        print(f"{'='*60}")
        
        results = {"status": "failed"}
        
        try:
            # Ensure we have a prepared protein to convert
            if not current_pdb_info.get("prepared_pdbqt"):
                print("❌ No prepared protein PDBQT file found. Skipping.")
                return results

            # Run PrankWeb prediction (CONSUME GENERATOR)
            generator = run_prankweb_prediction()
            final_output = None
            
            print("  ⚙️  Running PrankWeb (Selenium + OpenBabel)...")
            for output in generator:
                final_output = output
                if isinstance(output, tuple) and len(output) > 0:
                    status_update = output[0]
                    if isinstance(status_update, dict):
                        status_html = status_update.get('value', '')
                        if 'Converting' in status_html:
                            print("     -> Converting PDBQT to PDB...")
                        elif 'Uploading' in status_html:
                            print("     -> Uploading to PrankWeb (this may take time)...")
            
            if not final_output:
                results["error"] = "No output from PrankWeb"
                print(f"❌ PrankWeb returned no output")
                return results
            
            status_update = final_output[0]
            df_update = final_output[1] if len(final_output) > 1 else None
            status_html = status_update.get('value', '') if isinstance(status_update, dict) else str(status_update)
            
            if isinstance(df_update, dict):
                prankweb_df = df_update.get('value')
            else:
                prankweb_df = df_update
            
            if prankweb_df is not None and not prankweb_df.empty:
                # Save CSV to batch directory 
                csv_path = step_dir / "binding_sites.csv"
                prankweb_df.to_csv(csv_path, index=False)
                
                # CRITICAL: Update global config so Docking finds this CSV
                current_pdb_info["prankweb_csv"] = str(csv_path)
                
                results["status"] = "success"
                results["csv_file"] = str(csv_path)
                results["num_pockets"] = len(prankweb_df)
                results["status_message"] = status_html
                print(f"✅ Binding site prediction complete: {len(prankweb_df)} pockets found")
            else:
                results["error"] = "No binding sites predicted"
                results["status_message"] = status_html
                print(f"⚠️ No binding sites found")
            
        except Exception as e:
            results["error"] = str(e)
            results["traceback"] = traceback.format_exc()
            print(f"❌ Error in binding site prediction: {e}")
        
        return results
    
    def process_docking(self, step_dir):
        """Step 5: Molecular docking."""
        print(f"\n{'='*60}")
        print(f"STEP 5: Molecular Docking")
        print(f"{'='*60}")
        
        results = {"status": "failed"}
        
        try:
            # Check prerequisites
            if not current_pdb_info.get("prepared_pdbqt"):
                print("❌ Missing prepared PDBQT.")
                return results
            if not current_pdb_info.get("prankweb_csv"):
                print("❌ Missing PrankWeb CSV.")
                return results

            # Run docking (CONSUME GENERATOR)
            generator = run_molecular_docking()
            final_output = None
            
            print("  ⚙️  Running Vina Docking...")
            for output in generator:
                final_output = output
            
            if not final_output:
                results["error"] = "No output from docking"
                print(f"❌ Docking returned no output")
                return results
            
            # Extract results
            status_update = final_output[0]
            df_update = final_output[1] if len(final_output) > 1 else None
            pose_update = final_output[2] if len(final_output) > 2 else None
            
            status_html = status_update.get('value', '') if isinstance(status_update, dict) else str(status_update)
            
            if isinstance(df_update, dict):
                summary_df = df_update.get('value')
            else:
                summary_df = df_update
            
            if isinstance(pose_update, dict):
                # Try getting choices or value
                pose_options = pose_update.get('choices', pose_update.get('value', []))
            else:
                pose_options = pose_update if pose_update else []
            
            if summary_df is not None and not summary_df.empty:
                # 1. Save Summary CSV to batch folder
                csv_path = step_dir / "docking_summary.csv"
                summary_df.to_csv(csv_path, index=False)
                results["csv_file"] = str(csv_path)
                print(f"  ✓ Saved: docking_summary.csv")
                
                # 2. Copy the docked PDB/PDBQT files from DOCKING_RESULTS_DIR to batch folder
                # The docking module saves files in DOCKING_RESULTS_DIR/pdb and /pdbqt
                docking_pdb_dir = os.path.join(DOCKING_RESULTS_DIR, "pdb")
                docking_pdbqt_dir = os.path.join(DOCKING_RESULTS_DIR, "pdbqt")
                
                # Copy complex PDBs
                copied_pdbs = []
                if os.path.exists(docking_pdb_dir):
                    for pdb_file in glob.glob(os.path.join(docking_pdb_dir, "*complex.pdb")):
                        shutil.copy2(pdb_file, step_dir)
                        copied_pdbs.append(os.path.basename(pdb_file))
                
                # 3. Generate Screenshots for Top Poses
                pose_images = []
                # Ensure we handle pose options correctly (sometimes it's a list of strings)
                if isinstance(pose_options, list):
                    top_poses = pose_options[:5] # Top 5
                else:
                    top_poses = []

                for i, pose_name in enumerate(top_poses):
                    try:
                        # Get HTML viewer for this pose
                        # Note: This function looks in DOCKING_RESULTS_DIR, so we must run it before cleaning up
                        viewer_result = display_docked_structure(pose_name)
                        
                        viewer_html = None
                        if isinstance(viewer_result, dict):
                            viewer_html = viewer_result.get('value')
                        else:
                            viewer_html = viewer_result
                        
                        if viewer_html and not viewer_html.startswith('<div'):
                            safe_pose_name = "".join(c if c.isalnum() or c in ('-', '_') else '_' for c in pose_name[:40])
                            img_path = step_dir / f"pose_{i+1}_{safe_pose_name}.png"
                            self.save_3d_viewer_as_image(viewer_html, img_path, f"Pose {i+1}")
                            pose_images.append(str(img_path))
                            print(f"  ✓ Saved image: {img_path.name}")
                    except Exception as e:
                        print(f"  ⚠️ Could not save pose {i+1}: {e}")
                
                results["status"] = "success"
                results["pose_images"] = pose_images
                results["docked_files"] = copied_pdbs
                results["num_poses"] = len(summary_df)
                results["status_message"] = status_html
                print(f"✅ Docking complete: {len(summary_df)} poses found")
            else:
                results["error"] = "Docking failed or no results"
                results["status_message"] = status_html
                print(f"❌ Docking failed - no poses found")
            
        except Exception as e:
            results["error"] = str(e)
            results["traceback"] = traceback.format_exc()
            print(f"❌ Error in docking: {e}")
        
        return results
    
    def process_admet(self, step_dir):
        """Step 6: ADMET analysis."""
        print(f"\n{'='*60}")
        print(f"STEP 6: ADMET Analysis")
        print(f"{'='*60}")
        
        results = {"status": "failed"}
        
        try:
            # Run ADMET prediction
            result = run_admet_prediction()
            
            if result is not None:
                msg, df, csv_path = result
                
                if csv_path and os.path.exists(csv_path):
                    # Copy CSV
                    dest_csv = step_dir / "admet_results.csv"
                    shutil.copy2(csv_path, dest_csv)
                    results["csv_file"] = str(dest_csv)
                
                if df is not None and not df.empty:
                    # Also save as JSON for easier parsing
                    json_path = step_dir / "admet_results.json"
                    df.to_json(json_path, orient='records', indent=2)
                    results["json_file"] = str(json_path)
                
                results["status"] = "success"
                results["message"] = msg
                print(f"✅ ADMET analysis complete")
            else:
                results["error"] = "ADMET analysis failed - run docking first"
                print(f"⚠️ ADMET analysis skipped (no docking results)")
            
        except Exception as e:
            results["error"] = str(e)
            results["traceback"] = traceback.format_exc()
            print(f"❌ Error in ADMET analysis: {e}")
        
        return results
    
    def process_single_protein(self, protein_input):
        """Process a single protein through the complete pipeline."""
        print(f"\n{'#'*80}")
        print(f"# Processing: {protein_input}")
        print(f"{'#'*80}")
        
        # Create folder structure
        protein_dir, step_dirs = self.create_protein_folders(protein_input)
        
        # Initialize results tracking
        pipeline_results = {
            "protein_input": protein_input,
            "start_time": datetime.now().isoformat(),
            "steps": {}
        }
        
        # Step 1: Structure Search 
        step1_results = self.process_structure_search(
            protein_input, 
            step_dirs["01_structure_search"]
        )
        pipeline_results["steps"]["structure_search"] = step1_results
        
        if step1_results["status"] != "success":
            pipeline_results["pipeline_status"] = "failed_at_step_1"
            pipeline_results["end_time"] = datetime.now().isoformat()
            self.save_pipeline_summary(protein_dir, pipeline_results)
            return pipeline_results
        
        # Step 2: Ramachandran Analysis (Now with fallback)
        step2_results = self.process_ramachandran(
            step_dirs["02_ramachandran_analysis"]
        )
        pipeline_results["steps"]["ramachandran"] = step2_results
        
        # Step 3: Protein Preparation (Meeko)
        step3_results = self.process_protein_preparation(
            step_dirs["03_protein_preparation"]
        )
        pipeline_results["steps"]["protein_preparation"] = step3_results
        
        if step3_results["status"] != "success":
            pipeline_results["pipeline_status"] = "failed_at_step_3"
            pipeline_results["end_time"] = datetime.now().isoformat()
            self.save_pipeline_summary(protein_dir, pipeline_results)
            return pipeline_results
        
        # Step 4: Binding Site Prediction (PrankWeb)
        step4_results = self.process_binding_sites(
            step_dirs["04_binding_site_prediction"]
        )
        pipeline_results["steps"]["binding_sites"] = step4_results
        
        # Step 5: Molecular Docking (Vina) 
        step5_results = self.process_docking(
            step_dirs["05_molecular_docking"]
        )
        pipeline_results["steps"]["docking"] = step5_results
        
        # Step 6: ADMET Analysis
        step6_results = self.process_admet(
            step_dirs["06_admet_analysis"]
        )
        pipeline_results["steps"]["admet"] = step6_results
        
        # Finalize
        pipeline_results["pipeline_status"] = "completed"
        pipeline_results["end_time"] = datetime.now().isoformat()
        
        self.save_pipeline_summary(protein_dir, pipeline_results)
        
        print(f"\n{'='*60}")
        print(f"✅ Pipeline complete for {protein_input}")
        print(f"Results saved to: {protein_dir}")
        print(f"{'='*60}\n")
        
        return pipeline_results
    
    def save_pipeline_summary(self, protein_dir, results):
        """Save pipeline summary to JSON."""
        summary_path = protein_dir / "pipeline_summary.json"
        with open(summary_path, 'w') as f:
            json.dump(results, f, indent=2)
    
    def run_batch(self, protein_list):
        """Run pipeline for multiple proteins."""
        print(f"\n{'#'*80}")
        print(f"# BATCH PROCESSING: {len(protein_list)} proteins")
        print(f"# Output directory: {self.output_base_dir}")
        print(f"{'#'*80}\n")
        
        batch_results = {
            "batch_start_time": datetime.now().isoformat(),
            "total_proteins": len(protein_list),
            "proteins": {}
        }
        
        for i, protein_input in enumerate(protein_list, 1):
            print(f"\n[{i}/{len(protein_list)}] Processing: {protein_input}")
            
            try:
                result = self.process_single_protein(protein_input)
                batch_results["proteins"][protein_input] = result
            except Exception as e:
                print(f"❌ Fatal error processing {protein_input}: {e}")
                batch_results["proteins"][protein_input] = {
                    "status": "fatal_error",
                    "error": str(e),
                    "traceback": traceback.format_exc()
                }
        
        batch_results["batch_end_time"] = datetime.now().isoformat()
        
        # Save batch summary
        batch_summary_path = self.output_base_dir / f"batch_summary_{self.timestamp}.json"
        with open(batch_summary_path, 'w') as f:
            json.dump(batch_results, f, indent=2)
        
        # Generate CSV summary
        self.generate_batch_csv_summary(batch_results)
        
        print(f"\n{'='*80}")
        print(f"✅ BATCH PROCESSING COMPLETE")
        print(f"Summary saved to: {batch_summary_path}")
        print(f"{'='*80}\n")
        
        return batch_results
    
    def generate_batch_csv_summary(self, batch_results):
        """Generate CSV summary of batch results."""
        summary_data = []
        
        for protein_input, result in batch_results["proteins"].items():
            row = {
                "protein_input": protein_input,
                "pipeline_status": result.get("pipeline_status", "unknown")
            }
            
            if "steps" in result:
                steps = result["steps"]
                
                # Structure search
                if "structure_search" in steps:
                    ss = steps["structure_search"]
                    row["protein_name"] = ss.get("protein_name")
                    row["pdb_id"] = ss.get("pdb_id")
                    row["structure_search_status"] = ss.get("status")
                
                # Other steps
                row["ramachandran_status"] = steps.get("ramachandran", {}).get("status")
                row["preparation_status"] = steps.get("protein_preparation", {}).get("status")
                row["binding_sites_status"] = steps.get("binding_sites", {}).get("status")
                row["docking_status"] = steps.get("docking", {}).get("status")
                row["admet_status"] = steps.get("admet", {}).get("status")
                
                # Counts
                if "binding_sites" in steps:
                    row["num_binding_pockets"] = steps["binding_sites"].get("num_pockets")
                
                if "docking" in steps:
                    row["num_docking_poses"] = steps["docking"].get("num_poses")
            
            summary_data.append(row)
        
        df = pd.DataFrame(summary_data)
        csv_path = self.output_base_dir / f"batch_summary_{self.timestamp}.csv"
        df.to_csv(csv_path, index=False)
        print(f"CSV summary saved to: {csv_path}")


# Example usage
if __name__ == "__main__":
    # Define your protein list
    protein_list = [
        "COX-2",
        "MEK-1",
        "CDK4",
        "STAT3",
        "SMAD4",
        "MCL-1",
        "PARP1",
        "EGFR",
        "CDK6",
        "IL-6",
        "PLA2",
        "ACC",
        "FASN",
        "CPT1A"

    ]
    
    # Create batch processor
    batch_processor = ProteinPipelineBatch(output_base_dir="batch_results")
    
    # Run batch processing
    results = batch_processor.run_batch(protein_list)
    
    print("\n" + "="*80)
    print("BATCH PROCESSING SUMMARY")
    print("="*80)
    for protein, result in results["proteins"].items():
        status = result.get("pipeline_status", "unknown")
        print(f"{protein}: {status}")
    print("="*80)