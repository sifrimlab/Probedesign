import sys
from pyfaidx import Fasta

transcript_file = sys.argv[1]
gencode_file = sys.argv[2]
transcripts = open(transcript_file,'r').readlines()
transcripts = [x.rstrip() for x in transcripts]
faidx_input = [t for t in transcripts]

if "." in faidx_input[0]:
    genes = Fasta(gencode_file,split_char="|",duplicate_action='first',)
else:
    genes = Fasta(gencode_file,split_char="|",duplicate_action='first',key_function = lambda x: x.split('.')[0])


with open("faidx_input.txt", "w") as textfile:
    for element in faidx_input:
        textfile.write(element + "\n")

with open("probes_for_alligning.fasta", "w") as textfile:
    for transcript in (faidx_input):
        try:
            textfile.write(f'>{genes[transcript].name}\n{genes[transcript]}\n')
        except:
            pass

with open("complement_initial_probes.fasta", "w") as textfile:
    for transcript in (faidx_input):
        try:
            textfile.write(f'>{genes[transcript].name}\n{genes[transcript][:].complement}\n')
        except:
            pass
