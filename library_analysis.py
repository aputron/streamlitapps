# %%
import streamlit as st
import pandas as pd
import csv
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import shutil
from matplotlib import transforms


# %%
def wig_reader(file):
    """Reading and formatting WIG file
    
    Parameters
    ---------
    file: str
        WIG file path

    Returns
    -------
    pandas.DataFrame
        Frequency of Nucleotides (order: A, C, G, T)

    """
    df = pd.read_csv(file, sep="\t", header=None, index_col=[0])
    df = df.drop(columns=[5, 6, 7])
    return df

#%%
def wt_ml(files, width, height):
    """Calculates and visualizes the variance in nucleotide frequency of a 
    mutant library using WT as a threshold to account for Nanopore sequencing and basecalling errors
    Requres WT filename to be inputted first
    
    Parameters 
    ----------
    files: List[str, str]
        Filenames of WT and ML WIG files respectively

    Returns
    -------
    seaborn.Figure
        Graph of nucleotide frequencies
    """
    # vafm = values apart from max
    
    # calculates WT nucleotide frequencies
    wt = wig_reader(files[0])
    vafm_wt = wt.apply(lambda row: row.drop(row.idxmax()), axis=1) 
    vafm_wt["mutations"] = vafm_wt.sum(axis=1)
    wt["variation"] = vafm_wt["mutations"] / wt.max(axis=1) * 100

    # calculates ML nucleotide frequencies
    ml = wig_reader(files[1])
    vafm_ml = ml.apply(lambda row: row.drop(row.idxmax()), axis=1)
    vafm_ml["mutations"] = vafm_ml.sum(axis=1)
    ml["variation"] = vafm_ml["mutations"] / ml.max(axis=1) * 100
    
    # calculates the difference in frequencies
    variance = []
    for i in range(len(ml)):
        ml_v = ml["variation"].iloc[i]
        wt_v = wt["variation"].iloc[i]
        variance.append(max(0, ml_v - wt_v))

    df_v = pd.DataFrame([variance])

    _fig, _ax = plt.subplots(figsize=(width, height))
    graph = sns.heatmap(df_v, cmap="magma", vmax=1, cbar=False, yticklabels=False)
    # graph = px.imshow(variance)

    return graph.get_figure()


def mutation_frequency(list_of_files, names):
    mutational_freq = pd.DataFrame()
    for files in range(len(list_of_files)):
        df = wig_reader(list_of_files[files])
        values_apart_from_max = df.apply(lambda row: row.drop(row.idxmax()), axis=1)
        values_apart_from_max["mutations"] = values_apart_from_max.sum(axis=1)

        sample_name = names[files]
        mutational_freq[sample_name] = (
            values_apart_from_max["mutations"]
            / values_apart_from_max["mutations"].max()
            * 100
        )

    _fig, _ax = plt.subplots(figsize=(35, 2))
    graph = sns.heatmap([mutation_frequency['High Background'][:563], mutation_frequency['Low Background'][:563],mutational_freq['High Signal'][:563],  mutational_freq['Low Signal'][:563]], cmap ='magma', vmax=90)
    graph = sns.heatmap(mutational_freq, cmap ='magma', vmax=90)
    #graph = px.imshow(mutational_freq, color_continuous_scale="delta", zmax=90)

    return graph


# %%
# Streamlit application
# Page title
st.markdown(
    """
    # Variation Visualization
    """
)

with st.expander("Mutant Library Variation", True):
    with st.form("ml-form"):
        file = st.file_uploader(
            "FILE UPLOADER: Input two WIG files, labelled as WT and ML",
            type=".wig",
            accept_multiple_files=True,
        )
        submitted = st.form_submit_button("Submit!")
        if submitted and len(file) == 2:
            st.write("UPLOADED!")
            os.mkdir("tempDir")
            files = []
            for uploaded_file in file:  # save the uploaded file remotely to make appropriate graphs easily
                with open(os.path.join("tempDir", uploaded_file.name), "wb") as f:
                    f.write(uploaded_file.getbuffer()[218:])
                    files.append(f"tempDir/{uploaded_file.name}")
            x = wt_ml(files, 10, 2)
            st.pyplot(x)
            shutil.rmtree("tempDir")
        elif submitted and len(file) != 2:
            error = ValueError("There should be 2 WIG files")
            st.exception(error)

with st.expander("Mutant Variation in sorted samples", True):
    with st.form("mv-form"):
        file = st.file_uploader(
            "FILE UPLOADER: Input your WIG files",
            type=".wig",
            accept_multiple_files=True,
        )
        submitted = st.form_submit_button("Submit!")
        if submitted and len(file) == 2:
            st.write("UPLOADED!")
            os.mkdir("tempDir")
            files = []
            labels = []
            for   uploaded_file in file:  # save the uploaded file remotely to make appropriate graphs easily
                with open(os.path.join("tempDir", uploaded_file.name), "wb") as f:
                    # Try to find "\n" in the buffer
                    f.write(uploaded_file.getbuffer()[230:])
                    files.append(f"tempDir/{uploaded_file.name}")
                    labels.append(uploaded_file.name)
            x = mutation_frequency(files, labels)
            st.plotly_chart(x)
            shutil.rmtree("tempDir")
