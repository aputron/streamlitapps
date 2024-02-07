#!/usr/bin/python3

import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO

def ref_align(bamfile, sample, output):
    """Extract the aligned regions and save them as fasta files

    Parameters:
    ----------
    bamfile: str
        Input BAM file
    sample: str
        Name of directory containing query sequence 
    output: str
        Output file name
    """

    bamfile = pysam.AlignmentFile(bamfile, 'rb')
    fr = []

    # reads in bamfile
    for b_read in bamfile:
        # get trimmed read
        read = Seq(b_read.query_sequence[b_read.qstart - 30:b_read.qend + 30])

        # save sequence info
        f_read = SeqRecord(read, id=b_read.query_name, description=f'sample={sample}')

        # get quality scores for corresponding DNA sequences
        fr.append(f_read)

    # write into a FASTA file
    Bio.SeqIO.write(fr, f"{output}", "fasta")

if __name__ == "__main__":
    from argparse import ArgumentParser
    
    parser = ArgumentParser(description='Extract reads by read name from bam file')
    parser.add_argument('-b','--bam', help='Bam file', required=True)
    parser.add_argument('-s', '--sample', help='Sample Name', required=True)
    parser.add_argument('-o', '--output', help="Output file path", required=True)
    options = parser.parse_args()
    ref_align(options.bam, options.sample, options.output)