import sys

transcript_file = sys.argv[1]
transcripts = open(transcript_file,'r').readlines()
transcripts = [x.rstrip() for x in transcripts]
faidx_input = [t+": " for t in transcripts]

with open("faidx_input.txt", "w") as textfile:
    for element in faidx_input:
        textfile.write(element + "\n")
