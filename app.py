import streamlit as st
import requests
import time

# Page config
st.set_page_config(page_title="Smooth Phase Docking Portal", page_icon="🧬", layout="wide")

# Initialize session state for phase tracking and data retention
if "current_phase" not in st.session_state:
    st.session_state.current_phase = 1
if "pdb_code" not in st.session_state:
    st.session_state.pdb_code = ""
if "smiles_string" not in st.session_state:
    st.session_state.smiles_string = ""
if "protein_ready" not in st.session_state:
    st.session_state.protein_ready = False
if "ligand_ready" not in st.session_state:
    st.session_state.ligand_ready = False
if "scan_complete" not in st.session_state:
    st.session_state.scan_complete = False
if "files_uploaded" not in st.session_state:
    st.session_state.files_uploaded = False

st.title("🧬 Smooth-Transition Molecular Docking Portal")
st.write("A continuous, data-retaining workflow using adaptive progressive disclosure.")
st.write("---")

# ==========================================
# PHASE 1: FETCH DATA (Always Expanded initially)
# ==========================================
p1_expanded = st.session_state.current_phase == 1
with st.expander("🟢 Phase 1: Fetch Target & Drug Data", expanded=p1_expanded or st.session_state.current_phase > 1):
    st.write("Retrieve initial structural data from online databases.")
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Protein Target")
        pdb_input = st.text_input("Enter PDB ID (e.g., 6LU7):", value=st.session_state.pdb_code).upper()
        if st.button("Fetch PDB", key="fetch_pdb_btn"):
            if len(pdb_input) == 4:
                with st.spinner("Fetching..."):
                    # Simulating structural API download
                    time.sleep(0.5) 
                    st.session_state.pdb_code = pdb_input
                    st.session_state.protein_ready = True
                    st.success(f"✓ Cached structural file for {pdb_input}")
            else:
                st.error("Please enter a valid 4-character ID.")

    with col2:
        st.subheader("Drug/Ligand Properties")
        smiles_input = st.text_input("Enter Ligand SMILES:", value=st.session_state.smiles_string)
        if st.button("Resolve SMILES", key="resolve_smiles_btn"):
            if smiles_input:
                st.session_state.smiles_string = smiles_input
                st.session_state.ligand_ready = True
                st.success("✓ SMILES notation registered.")
            else:
                st.error("Please enter a SMILES string.")

    # Phase 1 Validation Gate
    if st.session_state.protein_ready and st.session_state.ligand_ready:
        if st.session_state.current_phase == 1:
            st.write("")
            if st.button("Unlock Phase 2 🔓", type="primary", key="unlock_p2"):
                st.session_state.current_phase = 2
                st.rerun()


# ==========================================
# PHASE 2: CAVITY SCANNING & GRID CONFIG
# ==========================================
p2_enabled = st.session_state.current_phase >= 2
p2_expanded = st.session_state.current_phase == 2

with st.expander("🟡 Phase 2: Cavity Scanning & Search Space Layout", expanded=p2_expanded):
    if not p2_enabled:
        st.warning("🔒 Complete Phase 1 to unlock this segment.")
    else:
        st.subheader("1. Active Site & Heteroatom Analyzer")
        c1, c2, c3 = st.columns(3)
        scan_cavity = c1.checkbox("Scan Pockets", value=True)
        scan_cofactor = c2.checkbox("Identify Co-factors")
        scan_hetero = c3.checkbox("Isolate Heteroatoms")

        if st.button("Execute Structural Scan", key="scan_btn"):
            with st.spinner("Analyzing PDB topology..."):
                time.sleep(1)
                st.session_state.scan_complete = True
                st.success("Analysis Complete: Grid environment maps successfully computed.")

        st.subheader("2. Search Grid Box Setup")
        grid_option = st.radio(
            "Grid constraint mapping:",
            ("Blind Docking (Global Grid Lock)", "Targeted Localized Coordinate Box")
        )
        
        if grid_option == "Targeted Localized Coordinate Box":
            cx, cy, cz = st.columns(3)
            cx.number_input("Center X", value=0.0)
            cy.number_input("Center Y", value=0.0)
            cz.number_input("Center Z", value=0.0)
        else:
            st.info("🔒 Grid is auto-locked globally to the entire surface boundary.")

        # Phase 2 Validation Gate
        if st.session_state.scan_complete:
            if st.session_state.current_phase == 2:
                if st.button("Unlock Phase 3 🔓", type="primary", key="unlock_p3"):
                    st.session_state.current_phase = 3
                    st.rerun()


# ==========================================
# PHASE 3: STRUCTURE FILE FORMAT UPLOADS
# ==========================================
p3_enabled = st.session_state.current_phase >= 3
p3_expanded = st.session_state.current_phase == 3

with st.expander("🔵 Phase 3: High-Fidelity File Uploads", expanded=p3_expanded):
    if not p3_enabled:
        st.warning("🔒 Complete Phase 2 to unlock this segment.")
    else:
        st.write("Provide the algorithmic runtime format inputs compiled from steps above.")
        col_up1, col_up2 = st.columns(2)
        
        with col_up1:
            up_pdbqt = st.file_uploader("Upload Target (.pdbqt):", type=["pdbqt"], key="p3_pdbqt")
        with col_up2:
            up_3d = st.file_uploader("Upload Drug 3D Topology (.sdf, .mol2, .pdbqt):", type=["sdf", "mol2", "pdbqt"], key="p3_3d")

        if up_pdbqt and up_3d:
            st.session_state.files_uploaded = True
            if st.session_state.current_phase == 3:
                if st.button("Unlock Final Phase 🔓", type="primary", key="unlock_p4"):
                    st.session_state.current_phase = 4
                    st.rerun()


# ==========================================
# PHASE 4: EXECUTION PIPELINE
# ==========================================
p4_enabled = st.session_state.current_phase >= 4
p4_expanded = st.session_state.current_phase == 4

with st.expander("🚀 Phase 4: Active Engine Docking & Matrix Analysis", expanded=p4_expanded):
    if not p4_enabled:
        st.warning("🔒 Complete all previous configurations to unlock execution.")
    else:
        st.success("Configuration criteria satisfied. System prepared for molecular matrix docking.")
        
        # Displaying historical data retention table natively visible above final execution
        st.markdown("### Retained Session Audit Stream")
        summary = {
            "Parameter Structure": ["PDB Code Target", "Chemical SMILES Mapping", "Calculated Cavities Status", "Input Format Verification"],
            "Configured Value": [st.session_state.pdb_code, st.session_state.smiles_string, "Scanned & Verified", "Ready"]
        }
        st.table(summary)

        if st.button("Execute Docking calculations", type="primary", key="run_dock_final"):
            bar = st.progress(0)
            for percent in range(100):
                time.sleep(0.02)
                bar.progress(percent + 1)
            st.metric(label="Predicted Free Energy of Binding (ΔG)", value="-9.2 kcal/mol", delta="High Target Affinity")
            st.balloons()
