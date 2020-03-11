import csv
import random
from  Bio import SeqIO
import Bio
from Bio.Seq import Seq

readout_file = "/home/ceyhun/Readout_probes_information.csv"
barcode_file = "/home/ceyhun/barcodes_merfish.csv"



#
def find_barcode(transcript):
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
    #Transcript names are added to the barcoding file to get corresponding bits and RScodes
    transcripts = ["ENST00000421812.3","ENST00000325192.8"]
    for transcript in transcripts:
        print(find_barcode(transcript))

results = {}
fw_primer = Seq("GGGCCACGTCCTCAATCGAC")
rv_primer = Seq("CCCTCGCCAAGGTTCGCTAG")
with open(readout_file, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        #print(row)
        results[row["Readout probe name"]] = row["Sequence"]


readout_ids = find_barcode(transcript)
outputfasta = open("/home/ceyhun/encoding_probes.readoutids.fa","w")
for record in SeqIO.parse("/home/ceyhun/firstprobedesign_seq30mer.20overlap..split","fasta"):
        s = random.sample(readout_ids,3)

        readout_seq = [results[i] for i in s]
        record.seq = fw_primer +"A" + Seq(readout_seq[0]) + "A" + record.seq + "A" + Seq(readout_seq[1]) + Seq(readout_seq[2]) + "A" + rv_primer
        record.id =  record.id +"_" + "_".join(s)
        print(record.id)
        record.description = ""
        #print(record.seq)
        SeqIO.write(record, outputfasta, "fasta")
print("Encoding probes written to /home/ceyhun/encoding_probes.readoutids.fa")
