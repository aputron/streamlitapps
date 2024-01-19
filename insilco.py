#%%
import streamlit as st
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio.SeqIO
import random
import os
import shutil

#%%
def random_mutagenesis(template_dna, count):
    dna=list(str(template_dna))
    for i in range(min(count), max(count)*3):
        mutation=random.randrange(len(template_dna))
        dna[mutation]=random.choice(list('ATGC'))
    dna="".join(dna)
    return dna

def dif(a, b):
    return [i for i in range(len(a)) if a[i] != b[i]]

# ref_seq: Seq object
# overall_lib_size: integer, default=100000...usually a large number
# aa_sub_count: the number between which the aa substitutions should be present: list of 2 numbers
# filename: filename of the subsequent fasta files

def mutant_library(ref_seq, filename, overall_lib_size=500000, aa_sub_count=[4, 5] ):
    # translate reference sequence
    ref_aa=ref_seq.translate()
    
    # make a mutation library with the specified mutations range
    mutated_seqs=[]
    mut_aa=[]
    print(overall_lib_size)
    for i in range(overall_lib_size):
        mutated_seq=Seq(random_mutagenesis(ref_seq, aa_sub_count))
        mutated_seqs.append(mutated_seq)
        mutated_aa=mutated_seq.translate()  
        mut_aa.append(mutated_aa)
        
    # select for mutations with the required number of amino acid substitutions
    candidates_dna = []
    candidates_aa=[]
    for i in range(len(mut_aa)):
        count = sum(1 for a, b in zip(ref_aa, mut_aa[i]) if a != b)
        if len(aa_sub_count) ==1:
            if count == aa_sub_count[0]:
                candidates_dna.append(mutated_seqs[i])
                candidates_aa.append(mut_aa[i])
        else:
            if aa_sub_count[0]<=count<=aa_sub_count[1]:
                candidates_dna.append(mutated_seqs[i])
                candidates_aa.append(mut_aa[i])
    
    # save mutants as fasta file
    fasta_dna = []
    fasta_aa = []
    for i in range(len(candidates_dna)):
        diff = dif(ref_aa, candidates_aa[i])
        description=[]
        for x in diff:
            des = f"{ref_aa[x]}{x + 1}{candidates_aa[i][x]}"
            description.append(des)
        description = "|".join(description)
        dna_read = SeqRecord(Seq(candidates_dna[i]), id = "Mutant Library", description=description)
        fasta_dna.append(dna_read)
        aa_read = SeqRecord(Seq(candidates_aa[i]), id = "Mutant Library", description=description)
        fasta_aa.append(aa_read)
    Bio.SeqIO.write(fasta_dna, f"{filename}_dna.fasta", "fasta")
    Bio.SeqIO.write(fasta_aa, f"{filename}_aa.fasta", "fasta")
    
    return candidates_dna, candidates_aa

# Streamlit application
# Page title
st.markdown(
    """
    # ***In-silico*** Mutant Library Generation
    """
)

with st.form("ml-form"):
    file = st.file_uploader("FILE UPLOADER: Input the FASTA file of the DNA of your target", type='.fasta')
    st.markdown("**OR**")
    seq = Seq(st.text_area("Input your target DNA sequence"))
    num = st.slider("Number of mutations per mutant", 1, 20,(1, 4))
    tot = st.number_input("Number of mutants in the library", max_value=1000000)
    st.write(f"A library of {tot} mutants with {num} mutations in each mutant will be generated")
    submitted = st.form_submit_button("Submit!")
    if submitted :
        st.write("UPLOADED!")
        os.mkdir("tempDir")
        if file != None:
            with open(os.path.join("tempDir", file.name),"wb") as f:
                    f.write(file.getbuffer())
            for record in Bio.SeqIO.parse(f"tempDir/{file.name}", "fasta"):
                sequence = record.seq
            st.write(sequence)
        elif seq != None:
            sequence = seq
            st.write(sequence)
        mut_count = []
        for i in num:
            mut_count.append(i)
        aa, dna = mutant_library(sequence, "mutant_library", overall_lib_size=tot, aa_sub_count=mut_count)
        #shutil.rmtree("tempDir")
        st.write(aa)

#%%
seq="TCATTTCGGGTTAACAACAGAGTAGTTAACCGGAACGAAACGGTAACCTTTACCTTCCGCACGGATGTGACCGATACCCGGGAACGGTAAATGGGAAGCCGCGATCAGGTAACCACCTTTCGCCGCGTCCGCGAACGCTTTTTTACGTTCAACCGCCGCAGATTTACCGTCGATGTCCAGCTGGTTGGTAACAGACGGGTCGTCGAACTGAACCGCCGCAACCAGGATCAGGTCACCCAGCAGCGCCAGTTTCTGACCCTGAGATTCAACAACGTAGGTGGTGTGACCCGGGGTGTGGCCATGGGAAGCCAGCGCTTTGATACCCGGAACCAGGTCGGTGTTACCAGAGAACGGTTTGAATTTACCCGCTTTAACGTACGGATTTAAGGAAGCCATCGCACCTTTGAAGAAACCTTTAGATTCGTCGTCCGGCGCTTTGTCCAGGTTGGTCTGAGACAGCCAGAAGTCCGCTTCTTTCTGGTCCGCACGAACAACCGCGTTCGGGAACGCCAGCTGTTCACCAACCATCAGACCACCAACGTGGTCCGGGTGCATGTGGGTGATGTAGATTTCGTCAACCTGTTCCGGCTGGTAACCCGCCGCTTTCAGGTTCGCCGCCAGACGACCCAGGGTCGGACCGAACAGACCCGCCGCACCGGTGTCAACCAGAACCAGTTTAGAACCGGTGTTAACCAGGTAACCGGTAACAGAGGTTTCCAGCGGCGCTTTCTGGAAAGATTTCGCCAGCGCAGACTGGGTTTTCGGCGCCGGCTGGTTCAGACGTTTGTCAACCGGCAGCGCAACGGTACCGTCAGACAGCGCGGTGATTTCGAAGTCACCCAGCAGCATACGGTAGTAACCCGGCGCAGAGGTACGAACCTGCGGCGCCGCCGCAGACGCGTGGGTAACGAAAACCATCTGCGCCGCGGTGCACAGACCCGCCAGCAGCAGAGAGGTTTTGGTCAGGTTACGCAT"
rev = Seq(seq).reverse_complement()
print(rev)
# %%
