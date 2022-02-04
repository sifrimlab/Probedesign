import os
import sys
import yaml
from Bio import SeqIO

# Taking only probes who both passed the filtering step and are uniquely alligned

config = yaml.load(open(sys.argv[1],'r'), Loader=yaml.FullLoader)
baseDir = sys.argv[2]
outputDir = sys.argv[3]

filtered_probes = os.path.join(outputDir,config["output"]["filtered_probes"])
uniquely_mapped_probes_probenames = os.path.join(outputDir,config["output"]["uniquely_mapped_probes_probenames"])
uniquely_mapped_probes_filtered_final = config["output"]["uniquely_mapped_probes_filtered_final"]

selected=[]
with open(os.path.join(outputDir,config["output"]["uniquely_mapped_probes_probenames"])) as t:
    for line in t:
        line = line.strip()
        selected.append(line)

recordslist=[]
with open(uniquely_mapped_probes_filtered_final,"w") as outfile:
    for record in SeqIO.parse(config["output"]["filtered_probes"], "fasta"):
        if record.id.split("_RS")[0] in selected:
            print(record.id)
            record.description = ""
            recordslist.append(record)

    SeqIO.write(recordslist,outfile, "fasta")
