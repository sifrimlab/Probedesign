import sys
import Bio
import argparse
from Bio import SeqIO
import csv
from itertools import permutations
from random import sample
import random
#
with open("/home/ceyhun/Readout_probes_information.csv",mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        line_count = 0
        csv_reader.fieldnames
        print("OG dictionary is :"+str(csv_reader))
        for row in csv_reader:
            if line_count == 0:
                print(f'column names are: {", ".join(row)}')
                line_count += 1
            print(f'bit {row["Bit"]} \t{row["Readout probe name"]} with sequence: {row["Sequence"]} and dye :{row["Dye"]}')
            line_count += 1
        list1= (row["Readout probe name"])
        print(random.sample(list1,3))
        print(f'Procssed {line_count -1} Readouts')
#



with open("/home/ceyhun/Readout_probes_information.csv",mode='r') as f:
    csv_list = [[val.strip() for val in r.split(",")] for r in f.readlines()]
(_,*header), *data = csv_list
csv_dict = {}
for row in data:
    key, *values = row
    csv_dict[key] = {key: value for key, value in zip(header,values)}
    print(row)
print(csv_dict.keys())
res = key, val = random.choice(list(csv_dict.items()))
print("random pair : "+str(res))
test = csv_dict.popitem(),csv_dict.popitem(),csv_dict.popitem()

print(sample(res,1))
print(" test random : " +str(test))
