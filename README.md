# Probedesign
MERFISH Probe design test

# first step of probe design take complement of transcripts of interest and output to fasta file. 
 
$ faidx --complement /home/ceyhun/gencode.v33.pc_transcripts.fixed.fa "ENST00000421812.3": "ENST00000325192.8": > Probe_complement_seq

$ pyfasta split -n1 -k30 -o20 Probe_complement_seq
# split the transcript in 30nt with 20nt overlap 
#output name is = 'Probe_complement_seq30mer.20overlap..split'

#then use output file for probedesign.py 
