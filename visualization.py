import base64
import json

def show_structure(pdb_text: str, pdb_id: str, protein_name: str) -> str:
    """
    Robust 3D visualization.
    - Uses vanilla JS (no jQuery dependency issues).
    - Highlights Ligands (Pink) vs Protein (Cartoon Spectrum).
    - Includes error reporting on the screen.
    """
    
    # Safely escape the PDB string for JavaScript
    pdb_json = json.dumps(pdb_text)
    
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
                background: rgba(255, 255, 255, 0.85); padding: 8px 12px;
                border-radius: 6px; font-family: sans-serif; font-size: 12px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.2); z-index: 100;
                pointer-events: none;
            }}
            .color-box {{ 
                display: inline-block; width: 10px; height: 10px; 
                margin-right: 6px; border-radius: 50%; 
            }}
        </style>
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    </head>
    <body>
        <div id="error-log"></div>
        
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
                    // 1. Check if library loaded
                    if (typeof $3Dmol === 'undefined') {{
                        throw new Error("3Dmol.js failed to load. Check internet connection.");
                    }}

                    // 2. Initialize Viewer
                    var element = document.getElementById('container');
                    var config = {{ backgroundColor: 'white' }};
                    var viewer = $3Dmol.createViewer(element, config);

                    // 3. Add Model
                    var pdbData = {pdb_json};
                    if (!pdbData || pdbData.length < 10) {{
                        throw new Error("PDB data is empty or invalid.");
                    }}
                    viewer.addModel(pdbData, "pdb");

                    // 4. Apply Styles
                    // A. Protein Style (Cartoon Rainbow)
                    // hetflag:false matches the standard amino acid chain
                    viewer.setStyle(
                        {{hetflag: false}}, 
                        {{cartoon: {{color: 'spectrum'}}}}
                    );

                    // B. Ligand Style (Pink Stick + Sphere)
                    // hetflag:true matches ligands, water, ions
                    // We combine stick and sphere to make it very visible
                    viewer.addStyle(
                        {{hetflag: true}}, 
                        {{
                            stick: {{colorscheme: 'magentaCarbon', radius: 0.25}},
                            sphere: {{scale: 0.3, color: 'magenta'}} 
                        }}
                    );

                    // 5. Render
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
    # Using srcdoc as a fallback if data URI fails in some strict browsers, though data URI is standard
    iframe = f'<iframe src="data:text/html;base64,{b64}" width="100%" height="500" frameborder="0" style="border: 1px solid #ccc; border-radius: 8px;"></iframe>'
    
    return iframe