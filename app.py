"""Gaussian Explorer: Streamlit application entry point."""
from pathlib import Path
import numpy as np
import streamlit as st
from scr.gaussian_io import parse_log
from scr import optimisation, vibrations, orbitals, calculation, structure_editor

st.set_page_config(page_title="Gaussian Explorer", page_icon="⚛️", layout="wide")

css_file = Path(__file__).parent / "style.css"
if css_file.is_file():
    st.markdown(f"<style>{css_file.read_text()}</style>", unsafe_allow_html=True)

st.title("GaussVue")
st.caption("Gaussian output analysis tool, and editing tool for Gaussian input files.")

with st.sidebar:
    logo = Path(__file__).parent / "images" / "BNNLab_v3.png"
    if logo.is_file():
        st.image(str(logo))
    st.header("Acknowledgements")
    st.write("This web tool was built to create, read, visualise and edit Gaussian input and output files.")
    st.header("Disclaimer")
    st.write("This software was developed by BNNLab, with all rights reserved. It is offered 'as is', without warranty of any kind, express or implied. The user assumes all risk for any malfunctions, errors, or damages resulting from the use of this software. The creator is not responsible for any direct or indirect loss arising from its use.")
    st.header("Input files")
    log_file = st.file_uploader("Gaussian output (.log / .out)", type=["log", "out"])
    fchk_file = st.file_uploader("Formatted checkpoint (.fchk / .fch)", type=["fchk", "fch"])
    cube_file = st.file_uploader("Orbital cube (.cube / .cub)", type=["cube", "cub"])
    st.divider()
    st.subheader("Viewer appearance")
    style = st.selectbox("Representation", ["Ball and stick", "Stick", "Sphere", "Line"])
    bg = st.color_picker("Background", "#AABDE1")
    spin = st.toggle("Rotate automatically", False)
    st.caption("3D rendering uses 3Dmol.js. Internet is needed to load its JavaScript library.")

if not any((log_file, fchk_file, cube_file)):
    st.info("Upload Gaussian output for analysis, or open Structure Editor to generate a molecule from SMILES.")

data = None
if log_file:
    try:
        data = parse_log(log_file.getvalue())
    except Exception as exc:
        st.error(f"Could not parse log: {exc}")

numbers = np.asarray(getattr(data, "atomnos", []), dtype=int) if data is not None else np.array([], dtype=int)
coords = np.asarray(getattr(data, "atomcoords", []), dtype=float) if data is not None else np.empty((0, 0, 3))
energies = np.asarray(getattr(data, "scfenergies", []), dtype=float) if data is not None else np.array([])
frequencies = np.asarray(getattr(data, "vibfreqs", []), dtype=float) if data is not None else np.array([])
displacements = np.asarray(getattr(data, "vibdisps", []), dtype=float) if data is not None else np.empty((0, 0, 3))

if data is not None:
    cols = st.columns(4)
    cols[0].metric("Atoms", len(numbers))
    cols[1].metric("Geometries", len(coords))
    cols[2].metric("Vibrational modes", len(frequencies))
    cols[3].metric("Final SCF energy", f"{energies[-1] / 27.211386245988:.7f} Eh" if len(energies) else "Not available")

tabs = st.tabs(["🧬 Structure & optimisation", "〰️ Vibrations", "☁️ Molecular orbitals", "📋 Calculation data", "✏️ Structure editor"])
with tabs[0]:
    optimisation.render(numbers, coords, energies, log_file, style, bg, spin)
with tabs[1]:
    vibrations.render(numbers, coords, frequencies, displacements, data, style, bg, spin)
with tabs[2]:
    orbitals.render(cube_file, fchk_file, numbers, coords, style, bg)
with tabs[3]:
    calculation.render(data)
with tabs[4]:
    structure_editor.render(numbers, coords, data, bg)

st.divider()
st.caption("Gaussian Explorer | Data processed locally by Streamlit/cclib; 3Dmol.js is loaded from a public CDN.")
