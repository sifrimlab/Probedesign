import sys
import Bio
import argparse
import timeit
import math
import re
import pandas
import gtfparse
import pyfaidx
from Bio import SeqIO
from pyfaidx import Fasta
from gtfparse import read_gtf

# Read in the gtf file

#df = read_gtf("/home/ceyhun/gencode.v33.chr_patch_hapl_scaff.basic.annotation.gff3")

# selecting only the genes / transcripts / exons / or only chr 1 in appropriate column
#df_genes = df[df["feature"] == "gene"]
#df_transcripts = df[df["feature"] == "transcript"]
#df_exon = df[df["feature"] == "exon"]
#df_genes_chr1 = df_genes[df_genes["seqname"] == "chr1"]


# https://pypi.org/project/pyfaidx/
# import fasta file for coding genes from gencode
genes = Fasta('/home/ceyhun/gencode.v33.pc_transcripts.fixed.fa')


##shows all possible keys for gene id's
genes.keys()
##Example key usage showing all NT in key use [:]'
genes['ENST00000361789.2'][:]


#making a dictionary from fasta file so every key is combined to value
record_dict = SeqIO.index('/home/ceyhun/gencode.v33.pc_transcripts.fixed.fa', "fasta")
print(record_dict["ENST00000361789.2"])

#how to call up multiple transcripts and get their sequence_
keys = ['ENST00000393001.1','ENST00000361789.2']

seqs = []
for key in keys:
    seqs.append(record_dict.get_raw(key).decode())

print(seqs)



# Taking fasta file and changing recods so that it only has the first description
# outputfasta = open("/home/ceyhun/gencode.v33.pc_transcripts.fixed.fa","w")
# seq_records = []
# for seq_record in SeqIO.parse("/home/ceyhun/gencode.v33.pc_transcripts.fa", "fasta"):
    # data = seq_record.id.split('|')
    # seq_record.id = data[0]
    # seq_record.name = data[0]
    # seq_record.description = data[0]
    # seq_records.append(seq_record)
#
