import streamlit as st
import os
import shutil
import re
from vina import Vina
import pandas as pd
from stmol import showmol
import py3Dmol

# --- PAGE SETUP ---
st.set_page_config(page_title="In-Silico AutoDock Vina Engine", page_icon="🧬", layout="wide")
st.title("🧬 Autonomous Molecular Docking Portal")
st.write("An interface powered by AutoDock Vina for localized or blind grid docking parameters.")

# --- UTILS & CORE ENGINE ---
def extract_vina_scores(output_pdbqt_path):
    """Parses Vina output PDBQT files to pull estimated binding affinities."""
    scores = []
    if not os.path.exists(output_pdbqt_path):
        return scores
        
    with open(output_pdbqt_path, "r") as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                # Matches energy value (kcal/mol)
                match = re.search(r"RESULT:\s*([-\d.]+)", line)
                if match:
                    scores.append(float(match.group(1)))
    return scores

def visualize_docking(receptor_path, ligand_path):
    """Renders receptor structure alongside predicted complex alignment."""
    with open(receptor_path, "r") as f:
        receptor_data = f.read()
    with open(ligand_path, "r") as f:
        ligand_data = f.read()

    view = py3Dmol.view(width=800, height=500)
    # Add receptor
    view.addModel(receptor_data, "pdbqt")
    view.setStyle({"cartoon": {"color": "spectrum"}})
    
    # Add ligand (Rendered as distinct stick topology)
    view.addModel(ligand_data, "pdbqt")
    view.setStyle({"model": 1}, {"stick": {"colorscheme": "cyanCarbon", "radius": 0.2}})
    
    view.zoomTo()
    showmol(view, height=500, width=800)

# --- SIDEBAR GRID BOX BINDING CONFIG ---
st.sidebar.header("🗺️ Grid Box Coordinate Configuration")
st.sidebar.markdown("Specify center search parameters and dimensions below:")

col_cx, col_cy, col_cz = st.sidebar.columns(3)
center_x = col_cx.number_input("Center X", value=0.0, step=1.0)
center_y = col_cy.number_input("Center Y", value=0.0, step=1.0)
center_z = col_cz.number_input("Center Z", value=0.0, step=1.0)

col_sx, col_sy, col_sz = st.sidebar.columns(3)
size_x = col_sx.number_input("Size X", value=20.0, min_value=5.0, step=1.0)
size_y = col_sy.number_input("Size Y", value=20.0, min_value=5.0, step=1.0)
size_z = col_sz.number_input("Size Z", value=20.0, min_value=5.0, step=1.0)

exhaustiveness = st.sidebar.slider("Exhaustiveness (Search Intensity)", min_value=1, max_value=16, value=4, step=1)
n_poses = st.sidebar.slider("Maximum Binding Poses to Return", min_value=1, max_value=10, value=5)

# --- WORKFLOW PHASES LAYOUT ---
st.header("Step 1: Input High-Fidelity Geometry Components")
c1, c2 = st.columns(2)

with c1:
    st.subheader("Target Structure")
    receptor_file = st.file_uploader("Upload Target Macromolecule (.pdbqt format required):", type=["pdbqt"])

with c2:
    st.subheader("Drug Candidate")
    ligand_file = st.file_uploader("Upload Rigid/Flexible Compound (.pdbqt format required):", type=["pdbqt"])

st.write("---")
st.header("Step 2: Simulation Dashboard Execution")

if receptor_file and ligand_file:
    # Set up safe, unique sandbox pathing per run
    work_dir = "docking_sandbox"
    os.makedirs(work_dir, exist_ok=True)
    
    rec_path = os.path.join(work_dir, "receptor.pdbqt")
    lig_path = os.path.join(work_dir, "ligand.pdbqt")
    out_path = os.path.join(work_dir, "output_poses.pdbqt")

    # Clear previous execution instances
    if st.button("🚀 Run Analytical Molecular Docking Engine", type="primary"):
        try:
            # Write streams locally to workspace
            with open(rec_path, "wb") as f:
                f.write(receptor_file.getbuffer())
            with open(lig_path, "wb") as f:
                f.write(ligand_file.getbuffer())
            
            with st.spinner("Initializing AutoDock Vina Engine force fields..."):
                # Initialize the Vina instance
                v = Vina(sf_name='vina', cpu=1)
                
                st.info("🔄 Binding receptor target maps...")
                v.set_receptor(rec_path)
                v.set_ligand_from_file(lig_path)
                
                st.info("📐 Locking box constraint geometry grids...")
                v.compute_vina_maps(
                    center=[center_x, center_y, center_z],
                    box_size=[size_x, size_y, size_z]
                )
                
                st.info("⚙️ Simulating Monte Carlo conformational searching algorithms...")
                # Run optimization docking computation
                v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
                
                # Output tracking lines
                v.write_poses(out_path, n_poses=n_poses, overwrite=True)
                st.success("✨ Molecular Docking Simulation Concluded Successfully!")
                
            # --- RESULTS PROCESSING ---
            st.write("---")
            st.header("Step 3: Post-Docking Analysis Metrics")
            
            scores = extract_vina_scores(out_path)
            if scores:
                res_df = pd.DataFrame({
                    "Conformational Mode": [f"Mode {i+1}" for i in range(len(scores))],
                    "Binding Affinity Score (kcal/mol)": scores
                })
                
                res_col, view_col = st.columns([2, 3])
                with res_col:
                    st.subheader("Scoring Index")
                    st.dataframe(res_df, use_container_width=True)
                    
                    # File export pipeline
                    with open(out_path, "rb") as file_out:
                        st.download_button(
                            label="📥 Download Docked Conformations (.PDBQT)",
                            data=file_out,
                            file_name="docked_output_poses.pdbqt",
                            mime="text/plain"
                        )
                
                with view_col:
                    st.subheader("3D Active-Site Verification Canvas")
                    visualize_docking(rec_path, out_path)
            else:
                st.error("Docking executed but no output structures were parsed correctly. Verify grid file boundaries.")
                
        except Exception as e:
            st.error(f"Execution Exception Context Triggered: {str(e)}")
            st.warning("Ensure coordinates accurately align inside target volume coordinates.")
else:
    st.info("💡 Complete Step 1 by uploading both structural files to activate the simulation console.")
