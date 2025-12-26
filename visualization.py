"""
3D structure visualization using 3Dmol.js
"""

import base64


def show_structure(pdb_text: str, pdb_id: str, protein_name: str) -> str:
    """Create 3D visualization HTML for PDB structure using base64 encoding."""
    pdb_escaped = pdb_text.replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$').replace('\r', '')
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
        <style>
            * {{
                margin: 0;
                padding: 0;
                box-sizing: border-box;
            }}
            body {{
                font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
                background: linear-gradient(135deg, #1e3c72 0%, #2a5298 100%);
                overflow: hidden;
            }}
            #viewer {{
                width: 100vw;
                height: 100vh;
                background: #0a0e27;
            }}
            .info-panel {{
                position: absolute;
                top: 20px;
                left: 20px;
                background: rgba(255, 255, 255, 0.95);
                padding: 16px 20px;
                border-radius: 12px;
                box-shadow: 0 8px 32px rgba(0,0,0,0.3);
                z-index: 100;
                max-width: 280px;
            }}
            .info-panel h2 {{
                margin: 0 0 2px 0;
                font-size: 20px;
                color: #1e3c72;
                font-weight: 700;
            }}
            .info-panel .subtitle {{
                font-size: 12px;
                color: #666;
                margin-bottom: 0;
                font-weight: 500;
            }}
            .controls {{
                position: absolute;
                bottom: 20px;
                left: 50%;
                transform: translateX(-50%);
                background: rgba(255, 255, 255, 0.95);
                padding: 12px 20px;
                border-radius: 12px;
                box-shadow: 0 8px 32px rgba(0,0,0,0.3);
                z-index: 100;
                display: flex;
                gap: 10px;
                align-items: center;
            }}
            .controls h3 {{
                margin: 0 12px 0 0;
                font-size: 14px;
                color: #1e3c72;
                font-weight: 700;
                white-space: nowrap;
            }}
            .controls button {{
                margin: 0;
                padding: 8px 16px;
                border: none;
                border-radius: 8px;
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                cursor: pointer;
                font-size: 12px;
                font-weight: 600;
                transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                box-shadow: 0 3px 10px rgba(102, 126, 234, 0.3);
                white-space: nowrap;
            }}
            .controls button:hover {{
                transform: translateY(-2px);
                box-shadow: 0 6px 20px rgba(102, 126, 234, 0.5);
            }}
            .controls button:active {{
                transform: translateY(0);
            }}
            .loading {{
                position: absolute;
                top: 50%;
                left: 50%;
                transform: translate(-50%, -50%);
                color: white;
                font-size: 18px;
                font-weight: 600;
                z-index: 50;
            }}
        </style>
    </head>
    <body>
        <div id="viewer"></div>
        <div class="loading" id="loading">Loading structure...</div>
        <div class="info-panel">
            <h2>{pdb_id}</h2>
            <div class="subtitle">{protein_name}</div>
        </div>
        <div class="controls">
            <h3>Visualization Style</h3>
            <button onclick="setCartoon()">Cartoon</button>
            <button onclick="setStick()">Stick</button>
            <button onclick="setSphere()">Sphere</button>
            <button onclick="setLine()">Line</button>
        </div>
        <script>
            let viewer;
            const pdbData = `{pdb_escaped}`;
            
            window.onload = function() {{
                try {{
                    const element = document.getElementById('viewer');
                    viewer = $3Dmol.createViewer(element, {{
                        backgroundColor: '#0a0e27'
                    }});
                    
                    viewer.addModel(pdbData, "pdb");
                    viewer.setStyle({{}}, {{'cartoon': {{'color': 'spectrum'}}}});
                    viewer.zoomTo();
                    viewer.render();
                    
                    document.getElementById('loading').style.display = 'none';
                }} catch(e) {{
                    console.error('Error loading structure:', e);
                    document.getElementById('loading').textContent = 'Error loading structure';
                }}
            }};
            
            function setCartoon() {{
                viewer.setStyle({{}}, {{'cartoon': {{'color': 'spectrum'}}}});
                viewer.render();
            }}
            
            function setStick() {{
                viewer.setStyle({{}}, {{'stick': {{'colorscheme': 'Jmol'}}}});
                viewer.render();
            }}
            
            function setSphere() {{
                viewer.setStyle({{}}, {{'sphere': {{'colorscheme': 'Jmol'}}}});
                viewer.render();
            }}
            
            function setLine() {{
                viewer.setStyle({{}}, {{'line': {{'colorscheme': 'chainHetatm'}}}});
                viewer.render();
            }}
        </script>
    </body>
    </html>
    """
    
    b64 = base64.b64encode(html_content.encode()).decode()
    iframe = f'<iframe src="data:text/html;base64,{b64}" width="100%" height="600" frameborder="0" style="border-radius: 12px; box-shadow: 0 8px 32px rgba(0,0,0,0.1);"></iframe>'
    
    return iframe
