import streamlit as st
import pandas as pd
import requests
from io import StringIO
import time
from Bio.PDB import PDBParser, PDBIO, Select
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem import AllChem
import streamlit.components.v1 as components

# --- Safe Graphics Import for Headless Linux Server Deployment ---
try:
    from rdkit.Chem import Draw
    RDKIT_DRAW_AVAILABLE = True
except ImportError:
    RDKIT_DRAW_AVAILABLE = False

# --- Streamlit Page Configuration ---
st.set_page_config(
    page_title="In Silico Molecular Docking Studio",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Helper Functions ---
def fetch_pdb_file(pdb_id):
    """Fetches PDB file contents from the RCSB protein data bank."""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    return None

class NonHeteroSelect(Select):
    """Biopython selection class to filter out heteroatoms (HETATM records)."""
    def accept_residue(self, residue):
        return residue.get_id()[0] == " "

def parse_heteroatoms(pdb_text):
    """Parses heteroatoms and co-factors directly from raw PDB text lines."""
    hetero_data = []
    for line in pdb_text.splitlines():
        if line.startswith("HETATM"):
            res_name = line[17:20].strip()
            chain_id = line[21].strip()
            res_seq = line[22:26].strip()
            atom_name = line[12:16].strip()
            # Prevent duplicate residue tracking
            if not any(d['Residue'] == res_name and d['ID'] == res_seq for d in hetero_data):
                hetero_data.append({
                    "Residue": res_name,
                    "Chain": chain_id,
                    "ID": res_seq,
                    "Type": "Co-factor / Heteroatom" if res_name not in ["HOH", "WAT"] else "Water Molecule"
                })
    return pd.DataFrame(hetero_data)

def render_3d_viewer(pdb_str, unique_key, ligand_smiles=None, style="cartoon", height=400):
    """Generates an inline HTML/JS canvas containing py3Dmol for 3D visualization."""
    style_opts = f"{{ {style}: {{color: 'spectrum'}} }}"
    
    ligand_js = ""
    if ligand_smiles:
        mol = Chem.MolFromSmiles(ligand_smiles)
        if mol:
            mol = Chem.AddHs(mol)
            try:
                AllChem.EmbedMolecule(mol)
                mol_block = Chem.MolToMolBlock(mol)
                cleaned_block = mol_block.replace('\n', '\\n').replace('\r', '')
                ligand_js = f"""
                var ligand_mol = msv.addModel({cleaned_block}, "sdf");
                msv.setStyle({{model: ligand_mol}}, {{stick: {{colorscheme: 'cyanCarbon'}} }});
                """
            except Exception:
                pass # Fallback if 3D embedding fails for complex SMILES

    cleaned_pdb = pdb_str.replace('\n', '\\n').replace('\r', '')
    viewer_id = f"3d_viewer_{unique_key}"
    
    html_content = f"""
    <div id="{viewer_id}" style="height: {height}px; width: 100%; position: relative;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
        var element = document.getElementById('{viewer_id}');
        var msv = $3Dmol.createViewer(element, {{backgroundColor: '#111217'}});
        var protein_mol = msv.addModel({cleaned_pdb}, "pdb");
        msv.setStyle({{model: protein_mol}}, {style_opts});
        {ligand_js}
        msv.zoomTo();
        msv.render();
    </script>
    """
    components.html(html_content, height=height+10)

# --- App State Initialization ---
if 'pdb_text' not in st.session_state:
    st.session_state.pdb_text = None
if 'pure_protein' not in st.session_state:
    st.session_state.pure_protein = None
if 'pdb_id' not in st.session_state:
    st.session_state.pdb_id = None
if 'smiles' not in st.session_state:
    st.session_state.smiles = ""
if 'ligand_props' not in st.session_state:
    st.session_state.ligand_props = None

# --- Main Dashboard Header ---
st.title("🧬 Multi-Phase In Silico Docking Workspace")
st.markdown("Automate receptor sanitization, ligand feature analysis, and grid-targeted compound docking in a seamless unified pipeline.")
st.divider()

# =====================================================================
# PHASE 1: PROTEIN PREPARATION
# =====================================================================
st.header("📍 Phase 1: Receptor Preparation")
col1, col2 = st.columns([1, 2])

with col1:
    pdb_input = st.text_input("Enter 4-Character PDB ID:", max_chars=4, placeholder="e.g., 1IEP").strip()
    fetch_btn = st.button("Fetch and Prepare Structure", type="primary")
    
    if fetch_btn and pdb_input:
        with st.spinner("Retrieving coordinate records from RCSB..."):
            raw_text = fetch_pdb_file(pdb_input)
            if raw_text:
                st.session_state.pdb_id = pdb_input.upper()
                st.session_state.pdb_text = raw_text
                
                # Convert to "Pure Protein" (Strip HETATMs)
                parser = PDBParser(QUIET=True)
                pdb_fh = StringIO(raw_text)
                structure = parser.get_structure(st.session_state.pdb_id, pdb_fh)
                
                io = PDBIO()
                io.set_structure(structure)
                out_stream = StringIO()
                io.save(out_stream, NonHeteroSelect())
                st.session_state.pure_protein = out_stream.getvalue()
                st.success(f"Successfully processed {st.session_state.pdb_id}!")
            else:
                st.error("Failed to discover PDB ID. Double-check your code entry.")

if st.session_state.pdb_text:
    render_mode = st.selectbox("3D Render Style", ["cartoon", "sphere", "line"], key="p1_style")
    c1, c2 = st.columns(2)
    
    with
