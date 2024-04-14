import streamlit as st

st.set_page_config(
    page_title="Drug Discovery",
    page_icon="💊",
)

st.write("""
    <div style="display: flex; justify-content: center;">
        <img src="https://i.imgur.com/DbxHi1L.jpeg" alt="Centered Image">
    </div>
""", unsafe_allow_html=True)
new_title = '<p style="text-align: center;font-family:Segoe UI Black; color:#6D59B9; font-size: 58px;">Multitask learning for Drug Solubility and Ligand Affinity prediction from SMILES</p>'
st.markdown(new_title, unsafe_allow_html=True)


st.sidebar.page_link("main.py", label="🏠 Home")
st.sidebar.page_link("pages/input_main.py", label="🔠 Manual Input")
st.sidebar.page_link("pages/compare.py", label="🆚 Compare")
st.sidebar.page_link("pages/about.py", label="🔎 About")

st.write("""
    <div style="display: flex; justify-content: center;">
        <img src="https://i.imgur.com/N9XjzM2.jpeg" alt="Centered Image">
    </div>
""", unsafe_allow_html=True)