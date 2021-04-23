import glob

configfile: "config.yaml"

transcript_file = config["files"]["transcripts"]
transcripts = open(transcript_file,'r').readlines()
transcripts = [x.rstrip() for x in transcripts]
faidx_input = [t+": " for t in transcripts]
github_test = "github test"

complement_sequence = config["output"]["complement_sequence"]
split_complement_sequence = config["output"]["split_complement_sequence"]
encoding_probes = config["output"]["encoding_probes"]
filtered_probes = config["output"]["filtered_probes"]
probe_stats = config["output"]["probe_stats"]

gencode_file = config["files"]["gencode"]
readout_file = config["files"]["readout"]
barcode_file = config["files"]["barcode"]


#Execs
faidx_exec = config["exec"]["faidx"]
pyfasta_exec = config["exec"]["pyfasta"]
python_exec = config["exec"]["python"]

#Params
kmer_length = config["parameters"]["kmer_length"]
kmer_overlap = config["parameters"]["kmer_overlap"]

#Paths:
src_path = config["path"]["src"]

rule all:
    input:
        complement_sequence,
        split_complement_sequence,
        encoding_probes,
        filtered_probes,
        probe_stats

rule complement_sequence:
    input:
        transcript_file
    output:
        complement_sequence
    shell:
        '''
        {faidx_exec} --complement {gencode_file} {faidx_input} > {output}
        '''

rule split_sequence:
    input:
        complement_sequence
    output:
        split_complement_sequence
    params:
        tmp_file=complement_sequence.replace(".fasta","")
    shell:
        '''
        {pyfasta_exec} split -n1 -k{kmer_length} -o{kmer_overlap} {input}
        mv {params.tmp_file}.split.{kmer_length}mer.{kmer_overlap}overlap.fasta  {output}
        '''
rule design_probes:
    input:
        split_complement_sequence,
        readout_file,
        barcode_file
    output:
        encoding_probes
    shell:
        '''
        {python_exec} {src_path}/probedesign.py config.yaml
        '''

rule filter_probes:
    input:
        encoding_probes
    output:
        filtered_probes,
        probe_stats
    shell:
        '''
        {python_exec} {src_path}/filtering_probes.py config.yaml
        '''
