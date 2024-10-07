import streamlit as st
import FlowCal
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import shutil

def facs_graph(files, labels, colour, width, height, y_lim, x_lim, bins, x_axis):
    """Make a graph in landscape mode

    Parameters:
    ----------
    files: List[str]
        List of fcs file paths
    labels: List[str]
        List of Labels for each sample (str)
    colour: List[str]
        List of colours to be assigned for the respective sample
    width: int
        Width of the graph
    height: int
        Height of the graph
    y_lim: Tuple[int, int]
        Range of the y limits
    x_lim: Tuple[int, int]
        Range of the x limits
    x_axis: bool
        Show x-axis label
    
    Returns
    -------
    seaborn.Figure
        Graph of inputed FCS files
    """
    tol_data = []
    mean_all = []

    # read all files
    for file in files:
        data = FlowCal.io.FCSData(file)
        data = FlowCal.transform.to_rfi(data)
        mean = FlowCal.stats.mean(data, channels=['FITC-A'])
        fitca = np.array(data[:, ['FITC-A']])
        tol_data.append(fitca)
        mean_all.append(mean)

    # transform the shape from 3-dimensional to 2-dimensional
    data = np.reshape(tol_data, (len(tol_data), len(tol_data[0]))).T

    # save as a dataframe
    df = pd.DataFrame(data, columns=labels)

    # plot the graph using matplotlib and seaborn
    _fig, _ax = plt.subplots(figsize=(width, height))
    _ax.get_legend().remove()
    graph = sns.histplot(df, multiple="layer", log_scale=True, element="poly", palette=colour, bins=bins)
    graph.set_ylim(y_lim)
    graph.set_xlim(x_lim)

    # plot the mean lines
    for i in range(len(labels)):
        plt.axvline(mean_all[i][0], color=colour[i])

    # to label the x axis
    if x_axis:
        plt.xlabel("Relative Fluoresence Units (RFU)")
    else:
        # plt.tick_params(labelleft=False, left=False)
        graph.set(xlabel=None)
    
    return graph.get_figure()

def facs_graph_por(files, labels, colour, width, height, y_lim, x_lim, bins, y_axis=True):
    """Make a graph in portrait mode

    Parameters:
    ----------
    files: List[str]
        List of fcs file paths
    labels: List[str]
        List of Labels for each sample (str)
    colour: List[str]
        List of colours to be assigned for the respective sample
    width: int
        Width of the graph
    height: int
        Height of the graph
    y_lim: Tuple[int, int]
        Range of the y limits
    x_lim: Tuple[int, int]
        Range of the x limits
    bins: int
        Number of bins of the histplot
    y_axis: bool
        Show y-axis label
    
    Returns
    -------
    seaborn.Figure
        Graph of inputed FCS files
    """
        
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
            # label_f = os.path.splitext(file)
            data_f = FlowCal.io.FCSData(file)
            data_f = FlowCal.transform.to_rfi(data_f)
            mean_f = FlowCal.stats.mean(data_f, channels=['FITC-A'])
            data_f = data_f[:, ['FITC-A']]
            data_1 = np.append(data_1, data_f) # have all the samples in a single column (long data)
            mean.append(mean_f)
    # mean.append(mean_wt)
    
    df = pd.DataFrame(data_1, columns=["FITC-A"]) # convert numpy array to a Pandas dataframe to use in seaborn
    
    # create a label for each of the samples in the concatentated df.
    label = []
    for i in labels:
        for _ in range(size):
            label.append(i)
    
    df["Label"] = label
    df_f = df.replace(0, np.nan)
    df_f.dropna()
    print(df_f)
       
    # set figure size
    _fig, _ax = plt.subplots(figsize=(width, height))
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
                        alpha=0.2, legend=False, bins=bins)
    
    # set graph configurations
    graph.set_ylim(y_lim)
    graph.set_xlim(x_lim)
    # plt.xticks([500])
    graph.set(xlabel=None)
    
    # plot the population means on the graph
    for i in range(len(labels)):
        plt.axhline(mean[i][0], color=colour[i])
    # plt.axhline(mean_wt, color="#00FFC3")

    # to label the y axis
    if y_axis:
        plt.ylabel("Relative Fluoresence Units (RFU)")
    else:
        # plt.tick_params(labelleft=False, left=False)
        graph.set(ylabel=None)
    
    return graph.get_figure()

# to clear out the session
def clear_all():
    """Clears the streamlit session"""
    st.session_state = False
    return

#%%

if __name__ == "__main__":
    # Streamlit application
    # Page title
    st.markdown(
        """
        # FACS Visualization
        """
    )

    # Color picker for the graph themes
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
        colours = [c1, c2, c3, c4, c5] # default colours

    # form to submit FACS files and generate their respective graphs
    with st.form("my-form"):
        file = st.file_uploader("FILE UPLOADER", type='.fcs', accept_multiple_files=True)
        themes = ["Low", "Medium", "High", "Neutral", "Background"]
        user_theme = st.multiselect("Choose your theme for your graphs in order", themes, default=themes)
        orientation = st.checkbox("Portrait graph")
        submitted = st.form_submit_button("Submit!")
        if submitted and file is not None:
            st.write("UPLOADED!")
            os.mkdir("tempDir")
            user_colours = []
            for theme in user_theme:
                x = themes.index(theme)
                user_colours.append(colours[x])
            st.write("To change the graph dimensions, choose the width and height on the slider below, and click on 'Submit!'(^) once more")
            col1, col2 = st.columns([1,4]) # to fix the graph size when seeing it on the UI
            with col2:
                files = []
                labels = []
                for uploaded_file in file: # save the uploaded file remotely to make appropriate graphs easily
                    with open(os.path.join("tempDir",uploaded_file.name),"wb") as f:
                        f.write(uploaded_file.getbuffer())
                        files.append(f"tempDir/{uploaded_file.name}")
                        labels.append(uploaded_file.name)
                if orientation:
                    with col1: # have different dimensional defaults for different graph orientations
                        width = st.slider("Graph width", 2, 10)
                        height = st.slider("Graph height", 5, 20, 10)
                        y_lim = st.slider("Y axis limits", 1, 30000, (1, 3000))
                        x_lim = st.slider("X axis limits", 1, 10000, (1, 6000))
                        bins = st.slider("Bins for histogram", 70, 200, 100)
                        y_axis=st.checkbox("Y axis label", True)
                    x = facs_graph_por(files, labels, user_colours, width, height, y_lim, x_lim, bins, y_axis)
                else:
                    with col1: # have different dimensional defaults for different graph orientations
                        width = st.slider("Graph width", 10, 20)
                        height = st.slider("Graph height", 1, 20)
                        y_lim = st.slider("Y axis limits", 1, 10000, (1, 500))
                        x_lim = st.slider("X axis limits", 1, 10000, (1, 3000))
                        bins = st.slider("Bins for histogram", 70, 200, 100)
                        x_axis = st.checkbox("X axis label", True)
                    x = facs_graph(files, labels, user_colours, width, height, y_lim, x_lim, bins, x_axis)
                st.pyplot(x)
            shutil.rmtree("tempDir")  # remove the temporary files folder to avoid overlapping queries

    # clear the session     
    st.button("Clear all", on_click=clear_all)
