import sys
import Bio
import argparse
import timeit
import numpy as np
import csv
import yaml
import os, os.path
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from csv import DictWriter
config = yaml.load(open(sys.argv[1],'r'), Loader=yaml.FullLoader)
baseDir = sys.argv[2]
outputDir = sys.argv[3]

# Filtering probes based on predefined values, creating a probe stats csv for downstream freedom of changing filters to achieve desired probes.


# Creating a Results folder if it doesn't exist yet
if not os.path.exists("Results/"):
     os.mkdir("Results/")

# Reading in the required config files
correctlength_seq = []
GClowerbound = config["parameters"]["GC_lower_bound"]
GCupperbound = config["parameters"]["GC_upper_bound"]

templowerbound = config["parameters"]["Temp_lower_bound"]
tempupperbound = config["parameters"]["Temp_upper_bound"]
probe_length = config["parameters"]["probe_length"]

outputfasta = config["output"]["filtered_probes"]
uniquely_mapped_probes_probenames = open(os.path.join(outputDir,config["output"]["uniquely_mapped_probes_probenames"]))
filtered_names = uniquely_mapped_probes_probenames.read().splitlines()
probe_stats = config["output"]["probe_stats"]

# Creating the probe stats csv
filters= {}

with open(probe_stats,"w", newline='') as filter_file:
    fieldnames =["probe_id","uniquely_mapped","probelength","GC","meltingtemp", "sequence"]
    csv_writer=csv.DictWriter(filter_file,fieldnames=fieldnames)
    csv_writer.writeheader()
    with open(outputfasta,"w") as outfile:
        for record in SeqIO.parse(os.path.join(outputDir,config["output"]["encoding_probes"]),"fasta"):
                filters["uniquely_mapped"] = True if record.id.split("_RS")[0] in filtered_names else False
                filters["probelength"] = len(record.seq)
                filters["probe_id"] = (record.id)
                filters["GC"] = (GC(record.seq))
                filters["meltingtemp"] = (mt.Tm_NN(record.seq))
                filters["sequence"] = (record.seq)
                csv_writer.writerow(filters)
                #print(filters)

                # If the probes confirm to chosen filters below write them out as filtered probes
                if (    len(record.seq) == probe_length
                    and GC(record.seq) > GClowerbound
                    and GC(record.seq) < GCupperbound
                    and float(('%0.2f' % mt.Tm_NN(record.seq))) > templowerbound
                    and float(('%0.2f' % mt.Tm_NN(record.seq))) < tempupperbound
                    ):
                    correctlength_seq.append(record)
                    SeqIO.write(record, outfile, "fasta")
