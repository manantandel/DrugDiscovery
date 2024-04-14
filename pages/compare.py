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
import deepchem as dc
import pickle
import pandas as pd
import gzip

model_list = ["AdaBoost", "Decision Tree", "Gradient Boosting", "Huber", "KNN", "Random Forest", "SVR", "XGBoost"]
key_list = ["AdaBoost", "DT", "GBR", "Huber", "KNN", "RF", "SVR", "XGB"]
fingerprint_list = ['Circular Fingerprint', 'MACCS Fingerprint', 'RDKit Fingerprint']
fingerprint_key_list = ["Circular", "MACCS", "RDKit"]

def iupac_to_smiles(iupac):
    iupac = iupac.lower()
    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(iupac) + '/smiles'
        ans = urlopen(url).read().decode('utf8')
        return ans
    except:
        return 'Did not work'


def check_convert_smiles(check_data):
    try:
        mol = Chem.MolFromSmiles(check_data)
        if mol is None:
            smiles = iupac_to_smiles(check_data)
            if smiles == "Did not work" or smiles=="":
                return False, None
            else:
                return True, smiles
        return True, check_data
    except Exception as e:
        return False, None
    

def draw_2d_molecule(smiles_data):
    mol = Chem.MolFromSmiles(smiles_data)

    # plt.figure(figsize=(4, 4))
    fig, ax = plt.subplots(figsize=(4, 4))
    if mol is not None:
        Draw.MolToMPL(mol)
        plt.title(smiles_data)
    else:
        plt.title("Invalid SMILES")

    plt.axis('off')
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    return plt


def draw_3d_molecule(smi, style='stick'):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    mblock = Chem.MolToMolBlock(mol)

    view = py3Dmol.view(width=400, height=400)
    view.addModel(mblock, 'mol')
    view.setStyle({style:{}})
    view.zoomTo()

    showmol(view)


def smiles_to_arr(smiles, fingerprint):
    if fingerprint == 0:
        featurizer = dc.feat.CircularFingerprint()
    elif fingerprint == 1:
        featurizer = dc.feat.MACCSKeysFingerprint()
    elif fingerprint == 2:
        featurizer = dc.feat.RDKitDescriptors()

    arr_smiles = featurizer.featurize(smiles)
    return arr_smiles


def predict_values(smiles, fingerprint, key):
    smiles = smiles_to_arr(smiles, fingerprint)

    fingerprint_selected = fingerprint_key_list[fingerprint]
    model_selected = key_list[key]
    location = f"Models\{fingerprint_selected}\{model_selected}{fingerprint_selected}.pkl.gz"

    with gzip.open(location, 'rb') as f:
        model = pickle.load(f)
    
    result_set = model.predict(smiles)
    result_set = result_set[0]

    return result_set[0], result_set[1], result_set[2], result_set[3], result_set[4]


st.set_page_config(
    page_title="Drug Discovery: Compare Molecules",
    page_icon="🆚",
    layout="wide"
)


st.write("""
    <div style="display: flex; justify-content: center;">
        <img src="https://i.imgur.com/DbxHi1L.jpeg" alt="Centered Image">
    </div>
""", unsafe_allow_html=True)
new_title = '<p style="text-align: center;font-family:Segoe UI Black; color:#6D59B9; font-size: 90px;">Compare Molecules</p>'
st.markdown(new_title, unsafe_allow_html=True)

st.sidebar.page_link("main.py", label="🏠 Home")
st.sidebar.page_link("pages/input_main.py", label="🔠 Manual Input")
st.sidebar.page_link("pages/compare.py", label="🆚 Compare")
st.sidebar.page_link("pages/about.py", label="🔎 About")

st.header("Select 2 elements for comparison")

col1, col2 = st.columns(2, gap="large")
st.markdown(
    """<style>
    div[class*="row-widget stRadio"] > label > div[data-testid="stMarkdownContainer"] > p {
        font-size: 30px;
        font-weight: bold;
    }
        </style>
        """, unsafe_allow_html=True) 
with col1:
    c1 = st.container(border=True)
    c2 = st.container(border=True)
    c5 = st.container(border=True)


    with c1:
        radio1 = st.radio("Chemical Molecule 1", ["Saved", "New"])

        smiles1 = ""

        if radio1 == "Saved":
            st.markdown(
    """<style>
    div[class*="row-widget stSelectbox"] > label > div[data-testid="stMarkdownContainer"] > p {
        font-size: 30px;
        font-weight: bold;
    }
        </style>
        """, unsafe_allow_html=True)
            file = open("analyzed_molecules.txt", "r")
            read_content = file.readlines()
            file.close()
            ele_list = [ele for ele in read_content if not ele.isspace()]
            dropbox1 = st.selectbox('Select an element', ele_list, key="d1")
            smiles1 = dropbox1
            
        elif radio1 == "New":
            st.markdown(
    """<style>
    div[class*="row-widget stSelectbox"] > label > div[data-testid="stMarkdownContainer"] > p {
        font-size: 30px;
        font-weight: bold;
    }
        </style>
        """, unsafe_allow_html=True)
            st.markdown("""
            <style>
            div[class*="stTextInput"] label p{
            font-size: 30px;
            font-weight: bold;
            }
                            </style>
            """, unsafe_allow_html=True)
            smiles_input_1 = st.text_input(label="Enter IUPAC name or SMILES string")

            check_input1, smiles_input1 = check_convert_smiles(smiles_input_1)

            if check_input1 == True:
                file = open("analyzed_molecules.txt", "a")
                file.write(smiles_input1+"\n")
                file.close()
                smiles1 = smiles_input1
            elif check_input1 == "":
                st.error("Input cannot be empty.")
            else:
                st.error("Unable to interpret input. Please provide a SMILES string or an IUPAC name.")

        fingerprint_select1 = st.selectbox("Select Fingerprint", fingerprint_list)
        model_select1 = st.selectbox("Select Model", model_list)

           
with col2:
    c3 = st.container(border=True)
    c4 = st.container(border=True)
    c6 = st.container(border=True)

    with c3:
                   
        radio2 = st.radio("Chemical Molecule 2", ["Saved", "New"])

        smiles2 = ""

        if radio2 == "Saved":
            st.markdown(
    """<style>
    div[class*="row-widget stSelectbox"] > label > div[data-testid="stMarkdownContainer"] > p {
        font-size: 30px;
        font-weight: bold;
    }
        </style>
        """, unsafe_allow_html=True)
            file = open("analyzed_molecules.txt", "r")
            read_content = file.readlines()
            file.close()
            ele_list = [ele for ele in read_content if not ele.isspace()]
            dropbox2 = st.selectbox('Select an element', ele_list, key="d2")
            smiles2 = dropbox2.rstrip("\n")

        elif radio2 == "New":
            st.markdown(
    """<style>
    div[class*="row-widget stSelectbox"] > label > div[data-testid="stMarkdownContainer"] > p {
        font-size: 30px;
        font-weight: bold;
    }
        </style>
        """, unsafe_allow_html=True)
            st.markdown("""
            <style>
            div[class*="stTextInput"] label p{
            font-size: 30px;
            font-weight: bold;
            }
                            </style>
            """, unsafe_allow_html=True)
            smiles_input_2 = st.text_input(label="Enter IUPAC name or SMILES string", key='t2')

            check_input2, smiles_input2 = check_convert_smiles(smiles_input_2)

            if check_input2 == True:
                file = open("analyzed_molecules.txt", "a")
                file.write(smiles_input2+"\n")
                file.close()
                smiles2 = smiles_input2
            elif check_input2 == "":
                st.error("Input cannot be empty.")
            else:
                st.error("Unable to interpret input. Please provide a SMILES string or an IUPAC name.")

        fingerprint_select2 = st.selectbox("Select Fingerprint", fingerprint_list, key="fs2")
        model_select2 = st.selectbox("Select Model", model_list, key="ms2")
    
subBut = st.button("Submit")
st.markdown("""
        <style>
                        .stButton {
            display: flex;
            justify-content: center;
        }
        .stButton > button p{
            justify-content: center;
            width: 400px;
            text-align: center;
            font-size: 25px;
            cursor: pointer;
        } </style>""", unsafe_allow_html=True)
if subBut:
    model_index1,fingerprint_index1 = model_list.index(model_select1), fingerprint_list.index(fingerprint_select1)
    aLogP1, BEI1, LE1, LLE1, SEI1 = predict_values(smiles1, fingerprint_index1, model_index1)

    model_index2,fingerprint_index2 = model_list.index(model_select2), fingerprint_list.index(fingerprint_select2)
    aLogP2, BEI2, LE2, LLE2, SEI2 = predict_values(smiles2, fingerprint_index2, model_index2)

    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    
    c2.subheader("General Properties", divider="violet")
    c2.metric("Molecular weight", value = round(Descriptors.MolWt(mol1),5), delta=round(Descriptors.MolWt(mol1) -Descriptors.MolWt(mol2),5))
    c2.metric("Number of Atoms",value = mol1.GetNumAtoms(), delta = (mol1.GetNumAtoms() - mol2.GetNumAtoms()))
    c2.metric("Number of Rings", value= Chem.rdMolDescriptors.CalcNumRings(mol1), delta=(Chem.rdMolDescriptors.CalcNumRings(mol1) - Chem.rdMolDescriptors.CalcNumRings(mol2)))
    c2.metric("Number of Rotatable Bonds", value= Descriptors.NumRotatableBonds(mol1), delta=(Descriptors.NumRotatableBonds(mol1)-Descriptors.NumRotatableBonds(mol2)))
    c2.metric("Number of Aromatic Rings", value= Descriptors.NumAromaticRings(mol1), delta=(Descriptors.NumAromaticRings(mol1)-Descriptors.NumAromaticRings(mol2)))
    c2.metric("Number of Aliphatic Rings", value= (Chem.rdMolDescriptors.CalcNumRings(mol1) - Descriptors.NumAromaticRings(mol1)), delta=((Chem.rdMolDescriptors.CalcNumRings(mol1) - Descriptors.NumAromaticRings(mol1))-(Chem.rdMolDescriptors.CalcNumRings(mol2) - Descriptors.NumAromaticRings(mol2))))
    c2.subheader("Predicted Properties", divider="violet")
    c2.metric("ALogP", value = round(aLogP1, 4), delta=round((aLogP1-aLogP2), 4))
    c2.metric("Ligand BEI", value=round(BEI1, 4), delta=round(BEI1 - BEI2, 4))
    c2.metric("Ligand LE", value=round(LE1, 4), delta=round(LE1 - LE2, 4))
    c2.metric("Ligand LLE", value=round(LLE1, 4), delta=round(LLE1 - LLE2, 4))
    c2.metric("Ligand SEI", value=round(SEI1, 4), delta=round(SEI1 - SEI2, 4))

    c2.subheader("2D Visualization", divider="violet")
    plot_mol1 = draw_2d_molecule(smiles1)
    c2.pyplot(plot_mol1)

    
    c4.subheader("General Properties", divider="violet")
    c4.metric("Molecular weight", value = round(Descriptors.MolWt(mol2),5), delta=round(Descriptors.MolWt(mol2) -Descriptors.MolWt(mol1),5))
    c4.metric("Number of Atoms",value = mol2.GetNumAtoms(), delta = (mol2.GetNumAtoms() - mol1.GetNumAtoms()))
    c4.metric("Number of Rings", value= Chem.rdMolDescriptors.CalcNumRings(mol2), delta=(Chem.rdMolDescriptors.CalcNumRings(mol2) - Chem.rdMolDescriptors.CalcNumRings(mol1)))
    c4.metric("Number of Rotatable Bonds", value= Descriptors.NumRotatableBonds(mol2), delta=(Descriptors.NumRotatableBonds(mol2)-Descriptors.NumRotatableBonds(mol1)))
    c4.metric("Number of Aromatic Rings", value= Descriptors.NumAromaticRings(mol2), delta=(Descriptors.NumAromaticRings(mol2)-Descriptors.NumAromaticRings(mol1)))
    c4.metric("Number of Aliphatic Rings", value= (Chem.rdMolDescriptors.CalcNumRings(mol2) - Descriptors.NumAromaticRings(mol2)), delta=((Chem.rdMolDescriptors.CalcNumRings(mol2) - Descriptors.NumAromaticRings(mol2))-(Chem.rdMolDescriptors.CalcNumRings(mol1) - Descriptors.NumAromaticRings(mol1))))
    c4.subheader("Predicted Properties", divider="violet")
    c4.metric("ALogP", value = round(aLogP2, 4), delta=round((aLogP2-aLogP1), 4))
    c4.metric("Ligand BEI", value=round(BEI2, 4), delta=round(BEI2 - BEI1, 4))
    c4.metric("Ligand LE", value=round(LE2, 4), delta=round(LE2 - LE1, 4))
    c4.metric("Ligand LLE", value=round(LLE2, 4), delta=round(LLE2 - LLE1, 4))
    c4.metric("Ligand SEI", value=round(SEI2, 4), delta=round(SEI2 - SEI1, 4))


    c4.subheader("2D Visualization", divider="violet")
    plot_mol2 = draw_2d_molecule(smiles2)
    c4.pyplot(plot_mol2)