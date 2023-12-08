#%%
import streamlit as st
import FlowCal
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff
import pandas as pd
import os
import altair as alt
import shutil

#%%%
#f = plt.figure(figsize=(10,5))
#a = f.add_subplot(1,1,1)
files = ["fcs/WT.fcs", "fcs/WT_phenol.fcs"]
x = []
for i in files:
    data = FlowCal.io.FCSData(i)
    data = FlowCal.transform.to_rfi(data)
    mean = FlowCal.stats.mean(data, channels=['FITC-A'])
    fitca = np.array(data[:, ['FITC-A']])
    x.append(fitca)

x = np.reshape(x, (len(x), len(x[0]))).T

group_labels = ['Group_1', 'Group_2']
x_df = pd.DataFrame(x, columns=group_labels)

graph = sns.histplot(x_df, multiple='layer', log_scale=True, element='poly')
graph.set_ylim(1.4, 2000)
graph.set_xlim(1.5, 1000)

#%%
# function to make FACS graphs
# files: [list] of fcs file paths (str)
# labels: [list] of Labels for each sample (str)
# colour: [list] of colours to be assigned for the respective sample
# y axis: [Boolean] referring to the presence of the y-axis

def facs_graph(files, labels, colour, y_axis=False, width=1):

    # reads first file
    data = FlowCal.io.FCSData(files[0])
    data = FlowCal.transform.to_rfi(data) # transforms arbitrary fluorescence units to Relative Fluoresence Units
    mean = FlowCal.stats.mean(data, channels=['FITC-A']) # calculates the population mean
    mean = [mean]
    data_1 = data[:, ['FITC-A']] # Extracts the channel with green fluorescence
    size = len(data_1)

    # if you have more fcs files to plot in the same graph
    if len(files) != 1:   
        for file in files[1:]:
            label_f = os.path.splitext(file)
            data_f = FlowCal.io.FCSData(file)
            data_f = FlowCal.transform.to_rfi(data_f)
            mean_f = FlowCal.stats.mean(data_f, channels=['FITC-A'])
            data_f = data_f[:, ['FITC-A']]
            data_1 = np.append(data_1, data_f) # have all the samples in a single column (long data)
            mean.append(mean_f)
    #mean.append(mean_wt)
    
    df = pd.DataFrame(data_1, columns=["FITC-A"]) # convert numpy array to a Pandas dataframe to use in seaborn
    
    # create a label for each of the samples in the concatentated df.
    label = []
    for i in labels:
        for x in range(size):
            label.append(i)
    
    df["Label"]=label
    df_f = df.replace(0, np.nan)
    df_f.dropna()
    print(df_f)
       
    #set figure size
    fig, ax = plt.subplots(figsize=(width, 10))
    # set graph styles
    sns.set_style("white")
    sns.set_style("ticks")
    sns.set_context("talk")

    # plot each species
    graph = sns.histplot(df_f, y="FITC-A", 
                        log_scale=True, 
                        element='poly',
                        palette=colour,
                        hue = "Label",
                        alpha=0.2, legend=False)
    
    # set graph configurations
    graph.set_ylim(1.4, 10000)
    graph.set_xlim(0, 10000)
    plt.xticks([500])
    graph.set(xlabel=None)
    
    # plot the population means on the graph
    for i in range(len(labels)):
        plt.axhline(mean[i][0], color=colour[i])
    #plt.axhline(mean_wt, color="#00FFC3")

    # to label the y axis
    if y_axis == True:
        plt.ylabel("Relative Fluoresence Units (RFU)")
    else:
        plt.tick_params(labelleft=False, left=False)
        graph.set(ylabel=None)
    
    return graph.get_figure()#, df_f

def clear_all():
    st.session_state=False
    return

#%%
st.markdown(
    """
    # FACS Visualization
    """
)

with st.expander("Choose the colours for each theme"):
    col1, col2, col3, col4, col5 = st.columns([1,1,1,1,1])
    with col1:
        c1 = st.color_picker("Low", "#DBBC0D")
    with col2:
        c2 = st.color_picker("Medium", "#B7DB00")
    with col3:
        c3 = st.color_picker("High", "#126E00")
    with col4:
        c4 = st.color_picker("Neutral", "#12528e")
    with col5:
        c5 = st.color_picker("Background", "#B0B0B0")
    colours = [c1, c2, c3, c4, c5]

with st.form("my-form"):
    file = st.file_uploader("FILE UPLOADER", type='.fcs', accept_multiple_files=True)
    themes = ["Low", "Medium", "High", "Neutral", "Background"]
    user_theme = st.multiselect("Choose your theme for your graphs in order", themes, default=themes)
    submitted = st.form_submit_button("Submit!")
    if submitted and file is not None:
        st.write("UPLOADED!")
        os.mkdir("tempDir")
        user_colours = []
        for theme in user_theme:
            x = themes.index(theme)
            user_colours.append(colours[x])
        col1, col2, col3, col4, col5= st.columns([1,1, 2, 1, 1])
        with col3:
            files = []
            labels = []
            for uploaded_file in file:
                with open(os.path.join("tempDir",uploaded_file.name),"wb") as f:
                    f.write(uploaded_file.getbuffer())
                    files.append(f"tempDir/{uploaded_file.name}")
                    labels.append(uploaded_file.name)
            x = facs_graph(files, labels, user_colours, width=2.5)
            st.pyplot(x)
        shutil.rmtree("tempDir")  
        
st.button("Clear all", on_click=clear_all)
