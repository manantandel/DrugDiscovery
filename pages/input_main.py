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
import pandas as pd
import pickle
import deepchem as dc
from streamlit_option_menu import option_menu
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

    plt.figure(figsize=(4, 4))
    if mol is not None:
        Draw.MolToMPL(mol)
        plt.title(smiles_data)
    else:
        plt.title("Invalid SMILES")

    plt.axis('off')
    return plt


def draw_3d_molecule(smi, style='stick'):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    mblock = Chem.MolToMolBlock(mol)

    view = py3Dmol.view(width=640, height=500)
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
    location = f".\Models\{fingerprint_selected}\{model_selected}{fingerprint_selected}.pkl.gz"

    with gzip.open(location, 'rb') as f:
        model = pickle.load(f)
    
    result_set = model.predict(smiles)
    result_set = result_set[0]

    return result_set[0], result_set[1], result_set[2], result_set[3], result_set[4]


def df_to_arr(dataframe, fingerprint):
    if fingerprint == 0:
        featurizer = dc.feat.CircularFingerprint()
    elif fingerprint == 1:
        featurizer = dc.feat.MACCSKeysFingerprint()
    elif fingerprint == 2:
        featurizer = dc.feat.RDKitDescriptors()

    df_featurized = featurizer.featurize(dataframe[0])
    return df_featurized


def predict_df_values(dataframe, model, fingerprint):
    fingerprint_selected = fingerprint_key_list[fingerprint]
    model_selected = key_list[model]
    location = f".\Models\{fingerprint_selected}\{model_selected}{fingerprint_selected}.pkl.gz"

    featurized_df = df_to_arr(dataframe, fingerprint)
    with gzip.open(location, 'rb') as f:
        model = pickle.load(f)
    
    result_set = model.predict(featurized_df)
    predict_df_temp = pd.DataFrame(result_set, columns=['ALogP', "BEI", "LE", "LLE", "SEI"])
    final_predicted_df = pd.concat([dataframe, predict_df_temp], axis=1)

    return final_predicted_df


def convert_df_to_csv(df):
  return df.to_csv().encode('utf-8')


st.set_page_config(
    page_title="Drug Discovery: Manual Input",
    page_icon="🔠",
    layout="wide"
)


st.write("""
    <div style="display: flex; justify-content: center;">
        <img src="https://i.imgur.com/Clh8Y75.jpeg" alt="Centered Image">
    </div>
""", unsafe_allow_html=True)
new_title = '<p style="text-align: center;font-family:Segoe UI Black; color:#6D59B9; font-size: 90px;">User Input</p>'
st.markdown(new_title, unsafe_allow_html=True)


st.sidebar.markdown(
    """<style>
    div[class*="row-widget stPageLink p"] > label > div[data-testid="stMarkdownContainer"] > p {
        font-size: 1.2rem;
        font-weight: bold;
    }
        </style>
        """, unsafe_allow_html=True)
st.sidebar.page_link("main.py", label="🏠 Home")
st.sidebar.page_link("pages/input_main.py", label="🔠 Manual Input")
st.sidebar.page_link("pages/compare.py", label="🆚 Compare")
st.sidebar.page_link("pages/about.py", label="🔎 About")

# tab1, tab2 = st.tabs(["Input Text", "Input from File"])

selected = option_menu(menu_title = None, options=["Text Input", "CSV Input"], 
                       icons=["input-cursor-text", "filetype-csv"], orientation="horizontal",
                       )

if selected == "Text Input":
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
    
    smiles_input = st.text_input(label="Enter IUPAC name or SMILES string")

    fingerprint_select = st.selectbox("Select Fingerprint", fingerprint_list)
    
    model_select = st.selectbox("Select Model", model_list)

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
    
    submit_button = st.button("Submit")

    check_input, smiles_input = check_convert_smiles(smiles_input)

    if submit_button:
        if check_input == True:
            file = open("analyzed_molecules.txt", "a")
            file.write(smiles_input+"\n")
            file.close()

            if smiles_input:
                mol = Chem.MolFromSmiles(smiles_input)

                with st.expander("Visualization", expanded=True):
                    st.subheader('Visualization')
                    sub_tab1, sub_tab2 = st.tabs(["2D", "3D"])

                    with sub_tab1:
                        col68, col69, col70 = st.columns([1, 2, 1])
                        with col69:
                            plot_mol = draw_2d_molecule(smiles_input)
                            st.pyplot(plot_mol)

                    with sub_tab2:
                        col1, col2, col3 = st.columns([1, 2, 1])
                        with col2:
                            draw_3d_molecule(smiles_input)

                model_index,fingerprint_index = model_list.index(model_select), fingerprint_list.index(fingerprint_select)

                aLogP, BEI, LE, LLE, SEI = predict_values(smiles_input, fingerprint_index, model_index)

                with st.expander("General Properties", expanded=True):
                    st.subheader('General properties')

                    col1, col2, col3 = st.columns(3)
                    col4, col5, col6 = st.columns(3)

                    col1.metric("Molecular weight", value = round(Descriptors.MolWt(mol),5))
                    col2.metric("Number of Atoms",value = mol.GetNumAtoms())
                    col3.metric("Number of Rings", value= Chem.rdMolDescriptors.CalcNumRings(mol))
                    col4.metric("Number of Rotatable Bonds", value= Descriptors.NumRotatableBonds(mol))
                    col5.metric("Number of Aromatic Rings", value= Descriptors.NumAromaticRings(mol))
                    col6.metric("Number of Aliphatic Rings", value= (Chem.rdMolDescriptors.CalcNumRings(mol) - Descriptors.NumAromaticRings(mol)))
                
                with st.expander("Predicted Values", expanded=True):
                    st.subheader("Predicted Values")
                    col7, col8,col9 = st.columns(3)
                    col10, col11, col12 = st.columns(3)

                    col7.metric("ALogP",value = round(aLogP,4))
                    col8.metric("Ligand Efficiency BEI",value = round(BEI,4))
                    col9.metric("Ligand Efficiency LE",value = round(LE,4))
                    col10.metric("Ligand Efficiency LLE",value = round(LLE,4))
                    col11.metric("Ligand Efficiency SEI",value = round(SEI,4))

                with st.expander("Final Result", expanded=True):
                    st.subheader("Final Result")
                    if 1.35 <= aLogP <= 1.8:
                        st.success("The molecule is a perfect candidate for drug discovery.")
                    elif -1 < aLogP < 5:
                        st.info("The molecule can be considered if there are less final candidates.")
                    else:
                        st.error("The molecule cannot be used for drug discovery.")
                    
            elif smiles_input == "":
                st.error("Input cannot be empty.")
        else:
            st.error("Unable to interpret input. Please provide a SMILES string or an IUPAC name.")

elif selected == "CSV Input":
    st.markdown(
    """<style>
    div[class*="row-widget stSelectbox"] > label > div[data-testid="stMarkdownContainer"] > p {
        font-size: 30px;
        font-weight: bold;
    }
        </style>
        """, unsafe_allow_html=True)
    
    st.markdown(
    """
    <style>
    div[class*="st-emotion-cache-1on073z e1b2p2ww0"] > label > div[data-testid="stMarkdownContainer"] > p {
        font-size: 30px;
        font-weight: bold;
    }
    </style>
    """,
    unsafe_allow_html=True
)
    
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
    
    st.markdown("""
        <style>
                        .stDownloadButton {
            display: flex;
            justify-content: center;
        }
        .stDownloadButton > button p{
            justify-content: center;
            width: 400px;
            text-align: center;
            font-size: 25px;
            cursor: pointer;
        } </style>""", unsafe_allow_html=True)
    
    with st.expander("What file to upload?", expanded=False):
        st.subheader("What file to upload?")
        st.write("The file should be a Comma Seperated Value(CSV) file. The file should contain only a column of SMILES data.")
    uploaded_file = st.file_uploader("Choose a CSV file", accept_multiple_files=False, type=['.csv'])
    fingerprint_select = st.selectbox("Select Fingerprint", fingerprint_list, key='s2')

    model_select = st.selectbox("Select Model", model_list, key='m2')

    submit_button = st.button("Submit", key='sub2')
    if submit_button:
        if uploaded_file is not None:
            model_index,fingerprint_index = model_list.index(model_select), fingerprint_list.index(fingerprint_select)
            df = pd.read_csv(uploaded_file ,header=None)
           
            with st.expander("Original Data", expanded=True):
                st.subheader("Original Data")
                c1,c2,c3 = st.columns([1,2,1])
                with c2:
                    c2.dataframe(df, column_config={"0":"SMILES"}, use_container_width=True)
            
            predicted_df = predict_df_values(df, model_index,fingerprint_index )
            with st.expander("Predicted Values", expanded=True):
                st.subheader("Predicted Values")
                c1,c2,c3 = st.columns([1,2,1])
                with c2:
                    c2.dataframe(predicted_df, column_config={"0":"SMILES"})

            file_download = predicted_df.to_csv("Predicted.csv")
            st.download_button(
                            label="Download data as CSV",
                            data=convert_df_to_csv(predicted_df),
                            file_name='Download.csv',
                            mime='text/csv',
                            )
            
