import os
import csv
import sys
import Bio
import yaml
import random
from  Bio import SeqIO
from Bio.Seq import Seq

# Designing the probes according to MERFISH paper https://doi.org/10.1038/s41598-018-22297-7

# function for finding the barcodes for a specified transcript in the barcode_merfish file
def find_barcode(transcript, barcode_file):
    readout_selected = []
    with open(barcode_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            #print(row)
            readout_selected = []
            if row["id"] == transcript:
                for (key, value) in row.items():
                    if value == "1":
                        readout_selected.append(key)
                return(readout_selected)
        return(False)

if __name__ == "__main__":

    # READ CONFIG
    config = yaml.load(open(sys.argv[1],'r'), Loader=yaml.FullLoader)
    baseDir = sys.argv[2]
    outputDir = sys.argv[3]
    split_complement_initial_probes = sys.argv[4]
    transcript_file = os.path.join(baseDir,config["files"]["transcripts"])
    barcode_file = os.path.join(baseDir,config["files"]["barcode"])
    readout_file = os.path.join(baseDir,config["files"]["readout"])
    fw_primer = Seq(config["primers"]["forward"])
    rv_primer = Seq(config["primers"]["reverse"])
    outputfasta = config["output"]["encoding_probes"]
#     split_complement_initial_probes = os.path.join(outputDir,config["output"]["split_complement_initial_probes"])

    results = {}
    with open(readout_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            #print(row)
            results[row["Readout probe name"]] = row["Sequence"]
    with open(outputfasta,"w") as outfile:
        for transcript in open(transcript_file, 'r'):
            transcript = transcript.rstrip()
            readout_ids = find_barcode(transcript, barcode_file)

        # For every transcript add 3 of the 4 readout_sequences, add the forward and reverse primer
        # Add 'A' spacers between the different parts
            for record in SeqIO.parse(split_complement_initial_probes,"fasta"):
                    if record.id.split("_")[0] == transcript:
                        random.seed(30)
                        s = random.sample(readout_ids,3)
                        readout_seq = [results[i] for i in s]
                        record.seq = fw_primer +"A" + Seq(readout_seq[0]) + "A" + record.seq + "A" + Seq(readout_seq[1]) + "A" + Seq(readout_seq[2]) + "A" + rv_primer
                        record.id =  record.id +"_" + "_".join(s)
                        # print(record.id)
                        record.description = ""
                        SeqIO.write(record, outfile, "fasta")
