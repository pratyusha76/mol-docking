import streamlit as st
import requests
import os

# Set up page configurations
st.set_page_config(page_title="Phase-Based AI Docking Portal", page_icon="🧬", layout="wide")

# Initialize session state variables for data retention across phases
if "phase" not in st.session_state:
    st.session_state.phase = 1
if "pdb_code" not in st.session_state:
    st.session_state.pdb_code = ""
if "smiles_string" not in st.session_state:
    st.session_state.smiles_string = ""
if "protein_file_path" not in st.session_state:
    st.session_state.protein_file_path = None
if "ligand_file_path" not in st.session_state:
    st.session_state.ligand_file_path = None
if "grid_type" not in st.session_state:
    st.session_state.grid_type = "Blind Docking"

# Helper function to move to the next phase
def next_phase():
    st.session_state.phase += 1

def prev_phase():
    st.session_state.phase -= 1

# Progress tracker banner at the top
st.title("🧬 Step-by-Step Molecular Docking Portal")
phases_desc = ["Phase 1: Fetch & Prepare Data", "Phase 2: Cavity & Grid Selection", "Phase 3: File Format Uploads", "Phase 4: Run Docking & Analysis"]
st.info(f"**Current Progress:** {phases_desc[st.session_state.phase - 1]}")
st.write("---")

# ==========================================
# PHASE 1: FETCH DATA (PDB ID & SMILES)
# ==========================================
if st.session_state.phase == 1:
    st.header("Phase 1: Fetch Target and Drug Data")
    st.write("Retrieve structural files using online databases.")

    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Protein Target Configuration")
        pdb_input = st.text_input("Enter 4-Character PDB ID (e.g., 6LU7):", value=st.session_state.pdb_code).upper()
        if st.button("Fetch Protein PDB"):
            if len(pdb_input) == 4:
                # Simulated/Actual Fetching Logic from RCSB PDB
                pdb_url = f"https://files.rcsb.org/download/{pdb_input}.pdb"
                response = requests.get(pdb_url)
                if response.status_code == 200:
                    st.session_state.pdb_code = pdb_input
                    # Save locally for next stages
                    st.session_state.protein_file_path = f"{pdb_input}.pdb"
                    with open(st.session_state.protein_file_path, "w") as f:
                        f.write(response.text)
                    st.success(f"Successfully fetched and cached Protein PDB: {pdb_input}")
                else:
                    st.error("Failed to fetch PDB file. Check the ID.")
            else:
                st.error("Please enter a valid 4-character PDB ID.")

    with col2:
        st.subheader("Drug/Ligand Configuration")
        smiles_input = st.text_input("Enter Ligand SMILES String:", value=st.session_state.smiles_string)
        if st.button("Resolve SMILES"):
            if smiles_input:
                st.session_state.smiles_string = smiles_input
                # In a real app, use RDKit to check validity or fetch properties from PubChem
                st.success("SMILES cached successfully for downstream 2D/3D arrangement.")
            else:
                st.error("Please enter a SMILES string.")

    st.write("---")
    # Condition to unlock Phase 2
    if st.session_state.protein_file_path and st.session_state.smiles_string:
        st.button("Proceed to Phase 2 ➡️", on_click=next_phase)
    else:
        st.warning("Please fetch both protein data and input a ligand SMILES to proceed.")


# ==========================================
# PHASE 2: SCAN CAVITY & LOCK GRID
# ==========================================
elif st.session_state.phase == 2:
    st.header("Phase 2: Cavity Scanning & Search Space Configuration")
    st.write(f"Working with Cached Protein: `{st.session_state.pdb_code}`")

    st.subheader("1. Elements to Identify & Scan")
    scan_cavity = st.checkbox("Scan Binding Cavities / Pockets", value=True)
    scan_cofactor = st.checkbox("Identify Co-factors (e.g., NAD, FAD, HEM)")
    scan_hetero = st.checkbox("Detect Heteroatoms & Water Molecules")

    if st.button("Run Scanning Computations"):
        with st.spinner("Analyzing structural properties..."):
            # Mock data processing representing parsing of HETATM lines in PDB
            st.info("Parsing structural metadata...")
            st.success("Analysis Complete: Identified 2 potential active pockets and 1 co-factor molecule.")

    st.subheader("2. Grid Box Setup (Lock Options)")
    grid_option = st.radio(
        "Select your Docking Grid configuration constraint:",
        ("Blind Docking (Lock Active Site to Whole Protein Grid)", "Targeted Box (Define Specific X, Y, Z coordinates)"),
        index=0 if st.session_state.grid_type == "Blind Docking" else 1
    )
    st.session_state.grid_type = grid_option

    if grid_option == "Targeted Box":
        col1, col2, col3 = st.columns(3)
        with col1: st.number_input("Center X", value=0.0)
        with col2: st.number_input("Center Y", value=0.0)
        with col3: st.number_input("Center Z", value=0.0)
    else:
        st.info("🔒 Grid is Locked globally. The entire protein space will act as the search canvas.")

    st.write("---")
    col_nav1, col_nav2 = st.columns([1,1])
    with col_nav1: st.button("⬅️ Back to Phase 1", on_click=prev_phase)
    with col_nav2: st.button("Proceed to Phase 3 ➡️", on_click=next_phase)


# ==========================================
# PHASE 3: FILE RE-UPLOAD & STRUCTURAL FILES
# ==========================================
elif st.session_state.phase == 3:
    st.header("Phase 3: Structural Formats & Uploads")
    st.write("Provide structural conversions required explicitly by automated structural docking algorithms.")

    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Protein Structural Formats")
        uploaded_pdbqt = st.file_uploader("Upload Target Protein in (.pdbqt) format:", type=["pdbqt"])
        if uploaded_pdbqt:
            st.session_state.protein_file_path = uploaded_pdbqt.name
            st.success(f"Target format tracked: {uploaded_pdbqt.name}")

    with col2:
        st.subheader("Drug Chemical Arrangements")
        uploaded_2d = st.file_uploader("Upload Drug Molecule 2D Representation (.sdf, .mol)", type=["sdf", "mol"])
        uploaded_3d = st.file_uploader("Upload Drug Molecule 3D Conformation (.pdb, .pdbqt, .mol2)", type=["pdb", "pdbqt", "mol2"])
        
        if uploaded_3d:
            st.session_state.ligand_file_path = uploaded_3d.name
            st.success("3D conformation tracked successfully.")

    st.write("---")
    col_nav1, col_nav2 = st.columns([1,1])
    with col_nav1: st.button("⬅️ Back to Phase 2", on_click=prev_phase)
    
    if uploaded_pdbqt and uploaded_3d:
        with col_nav2: st.button("Proceed to Phase 4 ➡️", on_click=next_phase)
    else:
        st.warning("Ensure both structural files (.pdbqt target and 3D compound file) are uploaded to move on.")


# ==========================================
# PHASE 4: EXECUTION & RETENTION REPORT
# ==========================================
elif st.session_state.phase == 4:
    st.header("Phase 4: Execution Pipeline & Analysis Data Retention")
    
    st.subheader("Summary of Structured Setup")
    summary_data = {
        "Parameter": ["Target PDB ID", "Ligand SMILES", "Grid Constraint Selection", "Target Implementation", "Drug Geometry File"],
        "Value/Status": [
            st.session_state.pdb_code if st.session_state.pdb_code else "Uploaded custom",
            st.session_state.smiles_string,
            st.session_state.grid_type,
            st.session_state.protein_file_path,
            st.session_state.ligand_file_path
        ]
    }
    st.table(summary_data)

    if st.button("🚀 Execute Molecular Docking Engine", type="primary"):
        progress_bar = st.progress(0)
        with st.spinner("Calculating configurations and energetic metrics..."):
            # Simulation of AutoDock Vina / AI Docking execution
            import time
            time.sleep(1.5)
            progress_bar.progress(40)
            time.sleep(1.5)
            progress_bar.progress(80)
            time.sleep(1)
            progress_bar.progress(100)
            
        st.success("✨ Calculation Matrix Rendered Successfully!")
        
        # Mocking Outputs metrics
        st.metric(label="Best Binding Affinity (kcal/mol)", value="-8.4", delta="Strong Binding Match")
        
        # Options to clear data session and start over cleanly
        if st.button("Clear Cache & Start New Run"):
            st.session_state.clear()
            st.rerun()

    st.write("---")
    st.button("⬅️ Back to Phase 3", on_click=prev_phase)
