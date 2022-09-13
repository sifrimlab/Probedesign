# Probedesign

## *Work in progress*

Takes as input transcript ids (ENST00000456328.2) and will result in a folder that contains a fasta file with the resulting filtered probes and probe stats including all the created probes with their values, in a table csv format.
The filtering requirements can be changed here aswell to achieve the desired probes. A ```Final_Probe_results.csv``` will give information about total created and filtered probes incuding how many probes per transcript are remaining.


Create conda environement using the provided probedesign.yml

Download reference transcripts __"wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.transcripts.fa.gz"__ and add them to the Inputs folder.

unzip the reference file __"gunzip gencode.v37.transcripts.fa.gz"__.

Generate latest reference genome index using STAR

In probedesign_processes.nf edit the params to match the location of the transcript files,genome index and specify an output folder.

To run the pipeline activate the conda environment and run: nextflow run probedesign_processes.nf
