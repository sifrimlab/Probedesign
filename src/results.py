import pandas as pd
from collections import Counter
import yaml
import sys
import csv
config = yaml.load(open(sys.argv[1],'r'), Loader=yaml.FullLoader)

initial_probes = 0
filtered_probes = 0
uniquely_alligned_probes = 0
filtered_and_alligned_probes = 0
alligned_transcripts=[]
filtered_transcripts=[]
filtered_and_alligned_transcripts=[]

with open(config["output"]["encoding_probes"]) as encoding:
    temp = encoding.read().splitlines()
    for line in temp:
        if "ENST" in line:
            initial_probes +=1
    # [initial_probes+=1 for line in temp if "ENST" in line]
    # initial_probes = sum([1 for line in temp if "ENST" in line])

with open(config["output"]["filtered_probes"]) as filtered1:
    temp = filtered1.read().splitlines()
    for line in temp:
        if "ENST" in line:
            filtered_probes +=1
            line = line.split("_RS")[0]
            filtered_transcripts.append(line)

    # [filtered_probes+=1 for line in temp if "ENST" in line]
    # [filtered_transcripts.append(line.split("_RS")[0]) for line in temp if "ENST" in line]

with open(config["output"]["uniquely_mapped_probes_probenames"]) as filtered2:
    temp = filtered2.read().splitlines()
    for line in temp:
        if "ENST" in line:
            line =">"+line.split("_RS")[0]
            uniquely_alligned_probes +=1
            alligned_transcripts.append(line)

    # [uniquely_alligned_probes+=1 for line in temp if "ENST" in line]
    # [alligned_transcripts.append(line.split("_RS")[0]) for line in temp if "ENST" in line]
for line in alligned_transcripts:
    if line in filtered_transcripts:
        line = line.split("_")[0]
        filtered_and_alligned_probes +=1
        filtered_and_alligned_transcripts.append(line)



duplicate_dict = Counter(filtered_and_alligned_transcripts)
df = pd.DataFrame(sorted(list(duplicate_dict.values()),reverse=True),index=duplicate_dict.keys())
df.to_csv("Results/Probes_per_gene.csv",header=False)

with open('Results/temp_results.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Initial probes",initial_probes])
    writer.writerow(["Uniquely alligned",uniquely_alligned_probes])
    writer.writerow(["Filtered probes",filtered_probes])
    writer.writerow(["Filtered and uniquely alligned",filtered_and_alligned_probes])
    writer.writerow(["Probes per transcript",""])

combined_csv = pd.concat([pd.read_csv("Results/temp_results.csv",header=None),pd.read_csv("Results/Probes_per_gene.csv",header=None)],axis=0,)
combined_csv.to_csv(config["output"]["Final_Probe_results"],index=None,header=None)
