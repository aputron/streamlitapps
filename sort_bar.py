import os
import pysam
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_list(filename):
    """Read the bamfile

    Parameters
    ----------
    filename: str
        Name of the bamfile

    Returns
    -------
    List[str]
        List of the read sequences
    """
    reads = []
    bamfile = pysam.AlignmentFile(filename, "rb")
    for b_read in bamfile:
        reads.append(b_read.query_sequence)
    return reads


default_barcodes_R = ["CCCTATGACA", "TAATGGCAAG", "AACAAGGCGT", "GTATGTAGAA", "TTCTATGGGG", "CCTCGCAACC", "TGGATGCTTA", "AGAGTGCGGC"]
default_barcodes_F = ["AATCCCACTAC", "TGAACTGAGCG", "TATCTGACCTT", "ATATGAGACG", "CGCTCATTAG", "TAATCTCGTC", "GCGCGATTTT", "AGAGCACTAG", "TGCCTTGATC", "CTACTCAGTC", "TCGTCTGACT", "GAACATACGG"]

def sort_reads(bamfile, folder_name, barcodes_R=None, barcodes_F=None):
    """Demultiplexes reads based on the combinations of forward and reverse barcode primers from a bamfile
    
    Parameters:
    ----------
    bamfile: str
        Name of the BAM file
    folder_name: str
        Path to directory containing query bamfile and output fasta files
    barcodes_R: List[str], optional
        Barcodes on the reverse primers (default is ["CCCTATGACA", "TAATGGCAAG", "AACAAGGCGT", "GTATGTAGAA", "TTCTATGGGG", "CCTCGCAACC", "TGGATGCTTA", "AGAGTGCGGC"])
    barcodes_F: List[str], optional
        Barcodes on the forward primers (default is ["AATCCCACTAC", "TGAACTGAGCG", "TATCTGACCTT", "ATATGAGACG", "CGCTCATTAG", "TAATCTCGTC", "GCGCGATTTT", "AGAGCACTAG", "TGCCTTGATC", "CTACTCAGTC", "TCGTCTGACT", "GAACATACGG"])
    """

    # defines the default barcodes
    if barcodes_R is None:
        barcodes_R = default_barcodes_R
    if barcodes_F is None:
        barcodes_F = default_barcodes_F
    
    # import the reads and their descriptions 
    reads = read_list(f'{folder_name}/{bamfile}')
    
    index, act_barcode = [], []
    for i in range(len(reads)):
        barcodes, bar = [], []
        for barcode in barcodes_R:
            rev = Seq(barcode).reverse_complement()
            if barcode.upper() in reads[i] or str(rev) in reads[i]:
                barcodes.append(barcodes_R.index(barcode) + 1)
                bar.append("R")
        for barcode in barcodes_F:
            rev = Seq(barcode).reverse_complement()
            if barcode.upper() in reads[i] or str(rev).upper() in reads[i]:   
                barcodes.append(barcodes_F.index(barcode) + 1)
                bar.append("F")
        bar = "_".join(bar)
        if bar == "R_F" and len(barcodes) == 2:
            act_bar = f'{barcodes[0]}_{barcodes[1]}'               
            with open(f'{folder_name}/{act_bar}.fasta', 'a') as f:
                f.write(f">{bar}\n")
                f.write(f'{reads[i]}\n')

sort_reads("barcode22-extracted-sorted.bam", "Georgii/georgii_231018_1/")

seqs = [] 
folder = "Georgii/georgii_231018_2"
subfolder = folder.split("/")[-1]
for file in os.listdir(f'{folder}/consensus'):
    if file.endswith("fasta"):
        filename = file.split(".")[0]
        for i in SeqIO.parse(f'{folder}/consensus/{file}', 'fasta'):
             dna_read = SeqRecord(i.seq, id=f'{filename}')
    seqs.append(dna_read)
SeqIO.write(seqs, f"{folder}/final/{subfolder}_consensus_96_seqs_dna.fasta", "fasta")