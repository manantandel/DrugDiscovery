import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors
from urllib.request import urlopen
from urllib.parse import quote
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
import py3Dmol
from stmol import showmol
from rdkit.Chem import AllChem


st.set_page_config(
    page_title="Drug Discovery",
    page_icon="💊",
)

st.write("""
    <div style="display: flex; justify-content: center;">
        <img src="https://i.imgur.com/GUUScJH.jpg" alt="Centered Image">
    </div>
""", unsafe_allow_html=True)

st.title("Drug Discovery")
st.sidebar.page_link("main.py", label="🏠 Home")
st.sidebar.page_link("pages/input_main.py", label="🔠 Manual Input")
st.sidebar.page_link("pages/compare.py", label="🆚 Compare")
st.sidebar.page_link("pages/about.py", label="🔎 About")
