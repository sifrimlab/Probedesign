import csv
import sys
import random
from  Bio import SeqIO
import Bio
from Bio.Seq import Seq
import yaml
config = yaml.load(open(sys.argv[1],'r'), Loader=yaml.FullLoader)

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
    # TODO Transcript names are added to the barcoding file to get corresponding bits and RScodes

    # READ CONFIG
    config = yaml.load(open(sys.argv[1],'r'), Loader=yaml.FullLoader)
    print(config)
    transcript_file = config["files"]["transcripts"]
    barcode_file = config["files"]["barcode"]
    readout_file = config["files"]["readout"]
    fw_primer = Seq(config["primers"]["forward"])
    rv_primer = Seq(config["primers"]["reverse"])
    outputfasta = open(config["output"]["encoding_probes"],"w")
    split_complement_initial_probes = config["output"]["split_complement_initial_probes"]

    results = {}
    with open(readout_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            #print(row)
            results[row["Readout probe name"]] = row["Sequence"]

    for transcript in open(transcript_file, 'r'):
        transcript = transcript.rstrip()
        readout_ids = find_barcode(transcript, barcode_file)


        for record in SeqIO.parse(split_complement_initial_probes,"fasta"):
                if record.id.split("_")[0] == transcript:
                    random.seed(30)
                    s = random.sample(readout_ids,3)
                    readout_seq = [results[i] for i in s]
                    record.seq = fw_primer +"A" + Seq(readout_seq[0]) + "A" + record.seq + "A" + Seq(readout_seq[1]) + "A" + Seq(readout_seq[2]) + "A" + rv_primer
                    record.id =  record.id +"_" + "_".join(s)
                    print(record.id)
                    record.description = ""
                    #print(record.seq)
                    SeqIO.write(record, outputfasta, "fasta")
