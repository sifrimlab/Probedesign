from Bio import SeqIO
import yaml
import sys

config = yaml.load(open(sys.argv[1],'r'), Loader=yaml.FullLoader)

filtered_probes = config["output"]["filtered_probes"]
uniquely_mapped_probes_probenames = config["output"]["uniquely_mapped_probes_probenames"]
uniquely_mapped_probes_selected_probes_fasta = open(config["output"]["uniquely_mapped_probes_selected_probes_fasta"],"w")

selected=[]
with open(config["output"]["uniquely_mapped_probes_probenames"]) as t:
    for line in t:
        line = line.strip()
        selected.append(line)

recordslist=[]
for record in SeqIO.parse(config["output"]["filtered_probes"], "fasta"):
    if record.id.split("_RS")[0] in selected:
        print(record.id)
        record.description = ""
        recordslist.append(record)

SeqIO.write(recordslist,uniquely_mapped_probes_selected_probes_fasta, "fasta")
