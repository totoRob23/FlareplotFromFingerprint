import pandas as pd
import json
import re
import random
import prolif as plf
from IPython.display import HTML
import colorsys

def generate_flareplot_from_fingerprint(fingerprint_path, threshold=0.5, width=500, height=500):
    """
    Generate an interactive flareplot from a fingerprint.pkl file
    
    Parameters:
    -----------
    fingerprint_path : str
        Path to the fingerprint.pkl file
    threshold : float, optional
        Threshold for filtering interactions (default: 0.5)
    width : int, optional
        Width of the flareplot (default: 500)
    height : int, optional
        Height of the flareplot (default: 500)
        
    Returns:
    --------
    IPython.display.HTML
        HTML object containing the interactive flareplot
    """
    # Load fingerprint and convert to dataframe
    def path_to_df(path):
        fp = plf.Fingerprint.from_pickle(path)
        return fp.to_dataframe()
    
    df = path_to_df(fingerprint_path)
    
    # Filter interactions based on threshold
    interaction_frequency = df.mean(axis=0)  
    filtered_columns = interaction_frequency[interaction_frequency >= threshold].index
    df_threshold = df[filtered_columns]
    
    # Function to generate a consistent color for a chain ID
    def get_chain_color(chain_id):
        # Use a hash of the chain ID to generate a consistent hue
        hue = (hash(chain_id) % 100) / 100.0  # Between 0 and 1
        saturation = 0.7
        value = 0.9
        
        # Convert HSV to RGB
        r, g, b = colorsys.hsv_to_rgb(hue, saturation, value)
        
        # Format as hex color
        return "#{:02x}{:02x}{:02x}".format(int(r*255), int(g*255), int(b*255))
    
    # Define ProLIF interaction colors
    interaction_colors = {
        "Hydrophobic": "#905000",
        "HBDonor": "#0000FF",
        "HBAcceptor": "#FF0000",
        "PiStacking": "#00FF00",
        "Anionic": "#A00000",
        "Cationic": "#0000A0",
        "CationPi": "#00A000",
        "PiCation": "#00A000",
        "VdWContact": "#505050",
        "EdgeToFace": "#00FF00",
        "FaceToFace": "#00FF00",
        "MetalAcceptor": "#8080FF",
        "MetalDonor": "#8080FF",
        "XBAcceptor": "#FF8080",
        "XBDonor": "#FF8080"
    }
    
    # Define amino acid three-letter to one-letter code mapping
    aa_code_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    # Function to convert ProLIF dataframe to enhanced JSON for flareplot
    def convert_prolif_df_to_enhanced_json(df):
        # Initialize the output structure with edges and tracks only
        output = {
            "edges": [],
            "tracks": [
                {
                    "trackLabel": "Boxes",
                    "trackProperties": []
                }
            ],
            "trees":  [
                {
                    "treeName": "Group-order",
                    "treeProperties": []
                }
            ]
        }
        
        # Get the column MultiIndex levels
        if isinstance(df.columns, pd.MultiIndex):
            ligands = df.columns.get_level_values('ligand').unique()
            proteins = df.columns.get_level_values('protein').unique()
            interactions = df.columns.get_level_values('interaction').unique()
        else:
            raise ValueError("Expected a DataFrame with MultiIndex columns")
        
        # Total number of frames
        total_frames = len(df.index)
        
        # Function to convert three-letter amino acid code to one-letter code
        def convert_aa_code(residue_string):
            # Extract the three-letter code and number from strings like "LYS123.A"
            match = re.match(r'([A-Z]{3})(\d+)\.([A-Z0-9]+)', residue_string)
            if match:
                three_letter = match.group(1)
                number = match.group(2)
                chain = match.group(3)
                # Convert to one-letter code if available
                one_letter = aa_code_map.get(three_letter.upper(), three_letter)
                return f"{chain}.{one_letter}{number}"
            return residue_string
        
        # Process each column (interaction)
        unique_nodes = set()
        chain_colors = {}  # Store colors for each chain
        
        # Process each column (interaction)
        for col_idx in range(len(df.columns)):
            ligand = df.columns[col_idx][0]
            protein = df.columns[col_idx][1]
            interaction_type = df.columns[col_idx][2]
            
            # Convert three-letter amino acid codes to one-letter codes
            name1 = convert_aa_code(ligand)
            name2 = convert_aa_code(protein)
            
            # Extract chain IDs
            chain1 = name1.split('.')[0]
            chain2 = name2.split('.')[0]
            
            # Generate colors for chains if not already done
            if chain1 not in chain_colors:
                chain_colors[chain1] = get_chain_color(chain1)
            if chain2 not in chain_colors:
                chain_colors[chain2] = get_chain_color(chain2)
            
            # Add nodes to tracking set
            unique_nodes.add(name1)
            unique_nodes.add(name2)
            
            # Count frames where this interaction exists
            frame_count = 0
            for frame_idx, frame_name in enumerate(df.index):
                if df.iloc[frame_idx, col_idx]:
                    frame_count += 1
            
            # Calculate weight as fraction of frames
            weight = frame_count / total_frames if total_frames > 0 else 0
            
            # Add to output if there are interactions
            if frame_count > 0:
                # Get color from ProLIF interaction colors or generate a fallback color
                edge_color = interaction_colors.get(interaction_type, "#{:06x}".format(hash(interaction_type) % 0xffffff))
                
                # Always use [0] for frames as requested
                output["edges"].append({
                    "name1": name1.split('.')[1], 
                    "name2": name2.split('.')[1],
                    "width": float(weight),
                    "frames": [0],
                    "color": edge_color 
                })
                
                output["tracks"][0]["trackProperties"].append({
                    "nodeName": name1.split('.')[1],
                    "color": chain_colors[chain1],
                    "size": float(weight)*2
                })
                
                output["tracks"][0]["trackProperties"].append({
                   "nodeName": name2.split('.')[1],
                   "color": chain_colors[chain2],
                   "size": float(weight)*2
                })
                
                output["trees"][0]["treeProperties"].append({"path": name1})
                output["trees"][0]["treeProperties"].append({"path": name2})
        
        return output
    
    # Convert the filtered dataframe to JSON
    json_data = convert_prolif_df_to_enhanced_json(df_threshold)
    
    # Create the HTML and JavaScript for the flareplot with SVG export functionality
    flareplot_html = f"""
    <div id="flare-container" style="width: {width}px; height: {height}px;"></div>
    <button id="export-svg">Export as SVG</button>
    <script>
        // Load D3.js
        require.config({{
            paths: {{
                d3: 'https://d3js.org/d3.v3.min'
            }}
        }});
        
        require(['d3'], function(d3) {{
            // Load flareplot script
            var script = document.createElement('script');
            script.src = 'https://cdn.rawgit.com/GPCRviz/FlarePlot/master/flareplot-main.js';
            script.onload = function() {{
                // Create the flareplot with data
                var jsonData = {json.dumps(json_data)};
                var plot = createFlareplot({width}, jsonData, "#flare-container");
                
                // Add export functionality
                d3.select("#export-svg").on("click", function() {{
                    // Get the SVG element
                    var svgElement = d3.select("#flare-container svg").node();
                    
                    // Get SVG source
                    var serializer = new XMLSerializer();
                    var source = serializer.serializeToString(svgElement);
                    
                    // Add namespaces
                    if(!source.match(/^<svg[^>]+xmlns="http\:\/\/www\.w3\.org\/2000\/svg"/)) {{
                        source = source.replace(/^<svg/, '<svg xmlns="http://www.w3.org/2000/svg"');
                    }}
                    if(!source.match(/^<svg[^>]+"http\:\/\/www\.w3\.org\/1999\/xlink"/)) {{
                        source = source.replace(/^<svg/, '<svg xmlns:xlink="http://www.w3.org/1999/xlink"');
                    }}
                    
                    // Add XML declaration
                    source = '<?xml version="1.0" standalone="no"?>\\r\\n' + source;
                    
                    // Convert SVG source to URI data scheme
                    var url = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(source);
                    
                    // Create download link
                    var link = document.createElement("a");
                    link.download = "flareplot.svg";
                    link.href = url;
                    document.body.appendChild(link);
                    link.click();
                    document.body.removeChild(link);
                }});
            }};
            document.head.appendChild(script);
        }});
    </script>
    """
    
    # Return the HTML object for display in the notebook
    return HTML(flareplot_html)

# Path to your fingerprint.pkl file
fingerprint_path = r"C:\Users\heiringr\Desktop\fingerprint.pkl"

# Generate and display the flareplot
flareplot = generate_flareplot_from_fingerprint(
    fingerprint_path, 
    threshold=0.5,  # Optional: adjust threshold for filtering interactions
    width=600,      # Optional: adjust width
    height=600      # Optional: adjust height
)

# Display the flareplot in the notebook
display(flareplot)

