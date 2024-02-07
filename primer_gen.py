#%%
import Bio #.SeqIO
import streamlit as st
import pandas as pd
import re

#%%

def read_mut_lib(filename):
    """Load the sequences and their descriptions

    Parameters:
    ----------
    filename: str
        FASTA file contianing the mutant library

    Returns
    -------
    Tuple[List[str], List[str]]
        List of Sequences and their corresponding descriptions
    """
    seq = []
    description = []    

    for record in Bio.SeqIO.parse(filename, "fasta"):
        seq.append(record.seq)
        description.append(record.description)
        
    return seq, description

# Function to extract the mutations to be made
# description: list of the descriptions of sequences in a FASTA file

def list_mut(description):
    """Extract the mutations to be made

    Parameters:
    ----------
    description: str
        Descriptions of sequences in a FASTA file      

    Returns
    -------
    List[List[str]]
        List of mutants metadata
    """
    des = []
    for i in description:
        string = i.split()
        mutant = string[1].split("|")
        des.append(mutant)
    return des

# Function to generate forward and reverse primers for a mutant library
def primer_db(filename):
    
    # Read the file and store sequences and descriptions
    seq, description = read_mut_lib(filename)
    
    # Lists that will be inputted into the resulting excel sheet
    ids=[] # mutant number
    for_primers = [] # forward primer
    rev_primers = [] # reverse primer
    mutant_f=[] # primer name (forward)
    mutant_r =[] # primer name (reverse) repeated for easily viewing on snapgene
    
    # Extract the mutations to be made
    des = list_mut(description)
    
    # identify primers and generate primers based on the mutations
    for i in range(len(des)):
        for x in des[i]:
            ids.append(i)
            mutant_f.append(x+"_F")
            mutant_r.append(x+"_R")
            
            # identify the aa number that is mutated
            numbers = re.findall(r'\d+', x)
            
            # identify the index of the dna that has to be changed
            base_change = int(numbers[0]) * 3
            
            # make a primer using 10 bases before and 9 bases after the mutable codon
            primer = seq[i][(base_change-12):(base_change+12)]
            
            # ensure there is no A at the 3' end
            if len(primer) != 0:
                if primer[-1]=="A":
                    count = base_change+12
                    while primer[-1]=="A":
                        primer+=seq[i][count]
                        count+=1
                elif primer[0]=="T":
                    count = (base_change-13)
                    while primer[0]=="T":
                        primer=seq[i][count]+primer
                        count-=1
                # make the reverse complement for the reverse primer 
                rev_primer = primer.reverse_complement()
            for_primers.append(str(primer))
            rev_primers.append(str(rev_primer))
            
    # Create a dataframe to easily make an excel sheet
    # Contains the mutant id, forward primer name, forward primer, reverse primer name, reverse primer
    df = pd.DataFrame({"Description":ids, "Primer_F":mutant_f, "Forward Primer":for_primers, "Primer_R":mutant_r, "Reverse Primer":rev_primers})
    
    return df

#%%
def valid_dna(seq):
    """Function to check if a given sequence has valid DNA nucleotides

    Parameters:
    ----------
    seq: str
        Query sequence

    Returns
    -------
    bool
        Does query sequence have the valid DNA nucleotides
    """
    dna_ref = set("ATGC")
    seq = set(seq.upper())
    return seq.issubset(dna_ref)

#%%

def generate_barcodes(forward, reverse, primer_name, forward_bc, reverse_bc):
    """Generate barcoded primers for multiplexed Nanopore sequencing

    Parameters:
    ----------
    forward: str
        5' to 3' region of the primer annealing to the forward strand
    reverse: str
        5' to 3' region of the primer annealing to the reverse strand
    primer_name: str
        General name of the primers; final name will include F or R and the number
    forward_bc: List[str]
        Barcodes to be attached to the forward primers; default primers used in the study
    reverse_bc: List[str]
        Barcodes to be attached to the reverse primers; default primers used in the study

    Returns
    -------
    pandas.DataFrame
        Barcoded primers with label
    """
    for_primers, for_primer_name, rev_primers, rev_primer_name = [], [], [], []

    # stich primers with barcodes
    for i in range(len(forward_bc)):
        for_primer_name.append(f'{primer_name}_F_{i + 1}')
        for_primers.append((forward_bc[i] + forward).upper())

    for i in range(len(reverse_bc)):
        rev_primer_name.append(f'{primer_name}_R_{i + 1}')
        rev_primers.append((reverse_bc[i] + reverse).upper())

    # make primer list
    primer_name = for_primer_name + rev_primer_name
    primers = {"Barcoded Primers": for_primers + rev_primers}

    return pd.DataFrame(primers, primer_name)

#%%

if __name__ == "__main__":
    # Streamlit application
    # Page title
    st.markdown(
        """
        # Primer Generation
        """
    )

    tab1, tab2 = st.tabs(["Mutant Primers", "Barcode Primers"])

    # to design primers for site-directed mutagenesis
    with tab1:
        with st.form("ml-form"):
            file = st.file_uploader(
            "FILE UPLOADER: Input one FASTA file containing your mutant candidates, with the respective mutations in the description",
            type=".fasta",
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
   
   
   
    # to generate barcoded primers
    with tab2:
        st.markdown("**Input your sequences (5' to  3'):**")
        col1, col2 = st.columns([1,1]) 
        with col1:
            forward = st.text_input("FORWARD") # input the forward annealing sequence
        with col2:
            reverse = st.text_input("REVERSE") # input the reverse annealing sequence

        with st.expander("Advanced"):
            col1, col2 = st.columns([1,1]) 
            with col1:
                primer_name = st.text_input("Primer names", value="Primer") # input general primer name
            with col2:
                filename = st.text_input("Barcode Filename", value="barcoded_primers") # input filname
            
            st.markdown("Inputting barcodes:")
            col1, col2 = st.columns([1,1]) 
            with col1:
                forward_bc_def = {"Forward Primers": [

                    "AATCCCACTAC", "TGAACTGAGCG", "TATCTGACCTT", 
                    "ATATGAGACG", "CGCTCATTAG", "TAATCTCGTC", 
                    "GCGCGATTTT", "AGAGCACTAG", "TGCCTTGATC", 
                    "CTACTCAGTC", "TCGTCTGACT", "GAACATACGG"
                    
                    ]} # default forward barcodes
                
                forward_bc = st.data_editor(forward_bc_def, num_rows="dynamic", use_container_width=True)
            
            with col2:
                reverse_bc_def = {"Reverse Primers": [

                    "CCCTATGACA", "TAATGGCAAG", "AACAAGGCGT", 
                    "GTATGTAGAA", "TTCTATGGGG", "CCTCGCAACC", 
                    "TGGATGCTTA", "AGAGTGCGGC"
                    
                    ]} # default reverse barcodes
                
                reverse_bc = st.data_editor(reverse_bc_def, num_rows="dynamic", use_container_width=True)
        
        submitted = st.button("Submit!")
        if submitted:
            # to check if the inputted barcodes are actually DNA bases
            f = "".join(forward_bc["Forward Primers"])
            r = "".join(reverse_bc["Reverse Primers"])

            # catch exceptions
            if forward == "" or reverse == "":
                st.error("Please input your forward/reverse sequences that you would like to stitch the barcodes to")
            elif not valid_dna(forward) or not valid_dna(reverse):
                st.error("Please input a valid DNA sequence")
            elif not valid_dna(f) or not valid_dna(r):
                st.error("Please input a valid DNA sequence as the barcode")
            else: # run script
                barcodes = generate_barcodes(forward, reverse, primer_name, forward_bc["Forward Primers"], reverse_bc["Reverse Primers"])
                st.table(barcodes)   # display table         
                st.download_button(label="Download primers", 
                                    data=barcodes.to_csv(), 
                                    file_name=f'{filename}.csv') # download excel file in an easy to order format
