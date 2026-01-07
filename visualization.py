import base64
import json
import os
from IPython.display import HTML

def show_structure(protein_text: str, ligand_text: str = None, pdb_id: str = "Docking Result", protein_name: str = "") -> str:
    """
    Robust 3D visualization supporting Overlay.
    
    Args:
        protein_text: The string content of the protein PDB (Rendered as Cartoon).
        ligand_text: (Optional) The string content of the ligand PDB (Rendered as Magenta Sticks).
        pdb_id: Title to display.
        protein_name: Subtitle to display.
    """
    
    # 1. Safely serialize the protein text
    prot_json = json.dumps(protein_text)
    
    # 2. Safely serialize the ligand text (if it exists)
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
                    // 1. Check Library
                    if (typeof $3Dmol === 'undefined') {{
                        throw new Error("3Dmol.js failed to load. Check internet connection.");
                    }}

                    // 2. Initialize Viewer
                    var element = document.getElementById('container');
                    var config = {{ backgroundColor: 'white' }};
                    var viewer = $3Dmol.createViewer(element, config);

                    // --- MODEL 0: PROTEIN ---
                    var proteinData = {prot_json};
                    if (!proteinData) throw new Error("Protein data is empty.");
                    
                    viewer.addModel(proteinData, "pdb");
                    
                    // Style Model 0 (Protein)
                    // We specifically target model 0 to keep styles separate
                    viewer.setStyle(
                        {{model: 0}}, 
                        {{cartoon: {{color: 'spectrum'}}}}
                    );
                    
                    // Optional: Show existing hetatoms in protein file as faint sticks
                    viewer.addStyle(
                        {{model: 0, hetflag: true}}, 
                        {{stick: {{radius: 0.1, color: 'lightgray'}}}}
                    );

                    // --- MODEL 1: LIGAND (OPTIONAL) ---
                    var ligandData = {lig_json};
                    
                    if (ligandData) {{
                        // Case A: Dual File Overlay
                        viewer.addModel(ligandData, "pdb");
                        
                        // Style Model 1 (Ligand)
                        // Target model 1 specifically
                        viewer.setStyle(
                            {{model: 1}}, 
                            {{stick: {{colorscheme: 'magentaCarbon', radius: 0.4}}}}
                        );
                    }} else {{
                        // Case B: Single File Fallback
                        // If no separate ligand file is provided, look for ligands INSIDE the protein file
                        viewer.addStyle(
                            {{resn: ["UNL", "LIG", "DRG", "UNK"]}}, 
                            {{stick: {{colorscheme: 'magentaCarbon', radius: 0.4}}}}
                        );
                    }}

                    // 4. Render
                    viewer.zoomTo();
                    viewer.render();

                }} catch(e) {{
                    logError(e.message);
                }}
            }};
        </script>
    </body>
    </html>
    """
    
    b64 = base64.b64encode(html_content.encode()).decode()
    iframe = f'<iframe src="data:text/html;base64,{b64}" width="100%" height="500" frameborder="0" style="border: 1px solid #ccc; border-radius: 8px;"></iframe>'
    
    return iframe

# --- HELPER FUNCTION TO READ FILES ---
def view_docking_files(protein_path, ligand_path):
    """Reads files and calls show_structure"""
    with open(protein_path, 'r') as f: p_text = f.read()
    with open(ligand_path, 'r') as f: l_text = f.read()
    return HTML(show_structure(p_text, l_text))