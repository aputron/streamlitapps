# %%
import streamlit as st
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio.SeqIO
import random
import sys
import os
from tqdm import tqdm
import shutil

#%%
def random_mutagenesis(template_dna, template_aa, mut_rng):
    """Generate random mutations across a template/target DNA

    Parameters
    ----------
    template_dna: Bio.SeqRecord.SeqRecord
        Template DNA to be mutated
    
    template_aa: Bio.SeqRecord.SeqRecord
        Template peptide sequence to be mutated

    mut_rng: Tuple[int, int]
        Range of mutations to be obtained
        
    Returns
    -------
    str
        Mutated DNA
    """
    template_aa = str(template_aa)
    dna = list(str(template_dna))
    no_mut = random.randrange(mut_rng[0], mut_rng[1] + 1)
    for i in range(no_mut):
        # Create only missense mutations
        condition = False
        while not condition:
            mutation = random.randrange(len(template_dna))
            mut_block_3 = mutation - mutation % 3
            dna[mutation] = random.choice("ATGC")
            aa_pos = mutation // 3 # floor division
            # translating the newly mutated sequence
            mutated_aa = Seq("".join(dna[mut_block_3:mut_block_3+3])).translate()
            if template_aa[aa_pos] != mutated_aa:
                condition = True # release the condition if the mutation is missense
            else:
                dna[mutation] = template_dna[mutation]
    return "".join(dna)

#%%
def dif(a, b):
    """Identifies the differences between two strings
    """
    return [i for i in range(len(a)) if a[i] != b[i]]

def mutant_library(ref_seq, filename, overall_lib_size=500000, aa_sub_count=(4, 5)):
    """Creates a mutant library of a specified size with specified number of mutations

    Parameters
    ----------
    ref_seq: Bio.Seq.Seq
        Reference Template DNA Sequence
    filename: str
        Name of file to save the resulting mutant library FASTA files
    overall_lib_size: int, optional
        Size of mutant library (default is 100000, usually a large number)
    aa_sub_count: Tuple[int, int]
        Range of the amino acid substitutions for each mutant

    Returns
    -------
    Tuple[List[Bio.Seq.Seq], List[Bio.Seq.Seq]]
        Mutant candidates' DNA and Amino Acid sequences
    """

    ref_aa = ref_seq.translate()

    # make a mutation library with the specified mutations range
    mutated_seqs, mut_aa = [], []
    for i in tqdm(range(overall_lib_size)):
        mutated_seq = Seq(random_mutagenesis(ref_seq, ref_aa, aa_sub_count))
        mutated_seqs.append(mutated_seq)
        mut_aa.append(mutated_seq.translate())

    # save mutants as fasta file
    fasta_dna = []
    fasta_aa = []
    # identify the mutation position and the aa substitution for the description
    for i in tqdm(range(len(mutated_seqs))):
        diff = dif(ref_aa, mut_aa[i])
        description = []
        for x in diff:
            des = f"{ref_aa[x]}{x + 1}{mut_aa[i][x]}"
            description.append(des)
        print(description)
        description = "|".join(description)
        dna_read = SeqRecord(
            Seq(mutated_seqs[i]), id="Mutant Library", description=description
        )
        fasta_dna.append(dna_read)
        aa_read = SeqRecord(
            Seq(mut_aa[i]), id="Mutant Library", description=description
        )
        fasta_aa.append(aa_read)
    Bio.SeqIO.write(fasta_dna, f"tempDir/{filename}_dna.fasta", "fasta")
    Bio.SeqIO.write(fasta_aa, f"tempDir/{filename}_aa.fasta", "fasta")

    return mutated_seqs, mut_aa

#%%
def valid_dna(seq):
    """Function to check if a given sequence has valid DNA nucleotides

    Parameters
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

if __name__ == "__main__":
    # Streamlit application
    # Page title
    st.markdown(
        """
        # ***In-silico*** Mutant Library Generation
        """
    )

    file = st.file_uploader(
        "FILE UPLOADER: Input the FASTA file of the DNA of your target", type=".fasta"
    )
    st.markdown("**OR**")
    seq = Seq(st.text_area("Input your target DNA sequence"))
    num = st.slider("Number of mutations per mutant", 1, 20, (1, 4))
    tot = st.number_input("Number of mutants in the library", max_value=1000000)
    st.write(
        f"A library of {tot} mutants with {num} mutations in each mutant will be generated"
    )
    submitted = st.button("Submit!")
    if submitted:
        shutil.rmtree("tempDir")
        st.write("UPLOADED!")
        os.mkdir("tempDir")
        if file is not None:
            with open(os.path.join("tempDir", file.name), "wb") as f:
                f.write(file.getbuffer())
            for record in Bio.SeqIO.parse(f"tempDir/{file.name}", "fasta"):
                sequence = record.seq
            st.write(sequence)
        elif seq is not None:
            sequence = seq
            st.write(sequence)
        if not valid_dna(sequence):
            st.error("That doesn't look like a DNA sequence. Please check your query once more")
        elif len(sequence) % 3 is not 0:
            st.error("The sequence cannot be translated. Please check the length of your input sequence")
        else:
            dna, aa = mutant_library(
                sequence, "mutant_library", overall_lib_size=tot, aa_sub_count=num
            )

            col1, col2 = st.columns([1,1]) 
            with col1:
                with open("tempDir/mutant_library_dna.fasta", "r") as dna_text:
                    st.download_button("DNA of Mutant candidates", data=dna_text, file_name='mutant_library_dna.fasta', mime='application/octet-stream')
                st.write(dna)
            with col2:
                with open("tempDir/mutant_library_aa.fasta", "r") as aa_text:
                    st.download_button("AA of Mutant candidates", data=aa_text, file_name='mutant_library_aa.fasta', mime='application/octet-stream')
                st.write(aa)
