# Probedesign
## TODO: Explain here how to get this from docker with conda environment

To run the pipeline run command ```snakemake --snakefile snakemake_probedesign```
For aditional snakemake command line specifications such as cores and threads check _https://snakemake.readthedocs.io/en/stable/executing/cli.html_.
In ```inputs/transcript_file.txt``` add the transcript id versions of interest -> Add these transcripts to ```Inputs/barcodes_merfish``` as well. 

Download reference transcripts __"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.transcripts.fa.gz"__ and add them to the Inputs folder
The filtering requirements and probe requirements can be changed  ```config.yaml``` as well as filenames and file outputs.
The snakemake pipeline can be adapted via ```snakefile_probedesign```.

The results folder contains a fasta file with the resulting filtered probes and probe stats including all the created probes with their values, in a table csv format.
The filtering requirements can be changed here aswell to achieve the desired probes. A ```Final_Probe_results.csv``` will give information about total created and filtered probes incuding how many probes per transcript are remaining.
