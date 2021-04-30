import pandas as pd
from collections import Counter
import yaml
import sys
import csv
config = yaml.load(open(sys.argv[1],'r'), Loader=yaml.FullLoader)
# Creating final report for designed probes

initial_probes = 0
filtered_probes = 0
uniquely_alligned_probes = 0
filtered_and_alligned_probes = 0
alligned_transcripts=[]
filtered_transcripts=[]
filtered_and_alligned_transcripts=[]

#Counting the probes/filtered probes/uniquely mapped probes

with open(config["output"]["encoding_probes"]) as encoding:
    temp = encoding.read().splitlines()
    for line in temp:
        if "ENST" in line:
            initial_probes +=1
  

with open(config["output"]["filtered_probes"]) as filtered1:
    temp = filtered1.read().splitlines()
    for line in temp:
        if "ENST" in line:
            filtered_probes +=1
            line = line.split("_RS")[0]
            filtered_transcripts.append(line)


with open(config["output"]["uniquely_mapped_probes_probenames"]) as filtered2:
    temp = filtered2.read().splitlines()
    for line in temp:
        if "ENST" in line:
            line =">"+line.split("_RS")[0]
            uniquely_alligned_probes +=1
            alligned_transcripts.append(line)

# Counting only probes who have both passed the filtering step and are uniquely alligned
for line in alligned_transcripts:
    if line in filtered_transcripts:
        line = line.split("_")[0]
        filtered_and_alligned_probes +=1
        filtered_and_alligned_transcripts.append(line)


# Counting how many times a probe exists for a specific transcript
duplicate_dict = Counter(filtered_and_alligned_transcripts)
df = pd.DataFrame(sorted(list(duplicate_dict.values()),reverse=True),index=duplicate_dict.keys())
df.to_csv("Results/Probes_per_gene.csv",header=False)

#Creating a temp results file with the count of the probes 
with open('Results/temp_results.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Initial probes",initial_probes])
    writer.writerow(["Uniquely alligned",uniquely_alligned_probes])
    writer.writerow(["Filtered probes",filtered_probes])
    writer.writerow(["Filtered and uniquely alligned",filtered_and_alligned_probes])
    writer.writerow(["Probes per transcript",""])

# Combining both csvs for a final probe results file including all the stats
combined_csv = pd.concat([pd.read_csv("Results/temp_results.csv",header=None),pd.read_csv("Results/Probes_per_gene.csv",header=None)],axis=0,)
combined_csv.to_csv(config["output"]["Final_Probe_results"],index=None,header=None)
