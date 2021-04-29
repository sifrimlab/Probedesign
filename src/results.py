import pandas as pd
from collections import Counter
import yaml
import sys
import csv
config = yaml.load(open(sys.argv[1],'r'), Loader=yaml.FullLoader)

initial_probes = 0
filtered_probes = 0
filtered_probes_unique_alligned = 0

with open(config["output"]["encoding_probes"]) as encoding:
    temp = encoding.read().splitlines()
    initial_probes = 0
    for line in temp:
        if "ENST" in line:
            initial_probes += 1

with open(config["output"]["filtered_probes"]) as filtered1:
    temp = filtered1.read().splitlines()
    for line in temp:
        if "ENST" in line:
            filtered_probes += 1

with open(config["output"]["filtered_probes_alligned_selected_sam"]) as filtered2:
    temp = filtered2.read().splitlines()
    for line in temp:
        if "ENST" in line:
            filtered_probes_unique_alligned += 1

with open(config["output"]["filtered_probes_alligned_selected_probenames"]) as f:
    temp = f.read().splitlines()
all_transcripts=[]
for line in temp:
    if "ENST" in line:
        line = line.split("_")[0]
        all_transcripts.append(line)

duplicate_dict = Counter(all_transcripts)
df = pd.DataFrame(sorted(list(duplicate_dict.values()),reverse=True),index=duplicate_dict.keys())
df.to_csv("Probes_per_gene.csv",header=False)

with open('results.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Initial probes",initial_probes])
    writer.writerow(["Filtered probes",filtered_probes])
    writer.writerow(["Filtered probes uniquely alligned",filtered_probes_unique_alligned])
    writer.writerow(["Probes per transcript",""])

combined_csv = pd.concat([pd.read_csv("results.csv",header=None),pd.read_csv("Probes_per_gene.csv",header=None)],axis=0,)
combined_csv.to_csv(config["output"]["Final_Probe_results"],index=None,header=None)
