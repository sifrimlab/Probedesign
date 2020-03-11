import sys
import Bio
import argparse
import timeit
import numpy as np
import csv
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from csv import DictWriter


correctlength_seq = []
GClowerbound = 43
GCupperbound = 63
outputfasta = open("/home/ceyhun/encoding_probes_filtered.fa","w")
filters= {}
with open("Filters.csv","w", newline='') as filter_file:
    fieldnames =["probe_id","probelength","GC","meltingtemp", "sequence"]
    csv_writer=csv.DictWriter(filter_file,fieldnames=fieldnames)
    csv_writer.writeheader()
    for record in SeqIO.parse("/home/ceyhun/encoding_probes.readoutids.fa","fasta"):
            filters["probelength"] = len(record.seq)
            filters["probe_id"] = (record.id)
            filters["GC"] = (GC(record.seq))
            filters["meltingtemp"] = (mt.Tm_NN(record.seq))
            filters["sequence"] = (record.seq)
            csv_writer.writerow(filters)
            #print(filters)

            if (    len(record.seq) == 134
                and GC(record.seq) > GClowerbound
                and GC(record.seq) < GCupperbound
                and float(('%0.2f' % mt.Tm_NN(record.seq))) > 66
                and float(('%0.2f' % mt.Tm_NN(record.seq))) < 76
                ):
                correctlength_seq.append(record)
                SeqIO.write(record, outputfasta, "fasta")
#print(filters)

print("Found %i correct sequences with a GC between 43%% and 63%% and melting temp between 66-76°C " % len(correctlength_seq))
print("Filters file written to Filters.csv")
print("Filtered encoding probes written to /home/ceyhun/encoding_probes_filtered")





# with open("/home/ceyhun/Probe_30nt20o.fa") as FASTA:
    # data = FASTA.read()
#
# print(data)
