# Probedesign

## *Work in progress*

Takes as input transcript ids (ENST00000456328.2) and will result in a folder that contains a fasta file with the resulting filtered probes and probe stats including all the created probes with their values, in a table csv format.
The filtering requirements can be changed here aswell to achieve the desired probes. A ```Final_Probe_results.csv``` will give information about total created and filtered probes incuding how many probes per transcript are remaining.


Create conda environement using the provided probedesign_nextflow.yml

Download reference transcripts __"wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.transcripts.fa.gz"__ and add them to the Inputs folder.

unzip the reference file in the inputs folder __"gunzip gencode.v37.transcripts.fa.gz"__.

Download the latest reference genome with accompanying GTF file in a folder :
mkdir genome_index 
cd genome_index 

__" wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"__ 
__" http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz"__
unzip these files with gunzip

Generate latest reference genome index using STAR.
__" STAR --runThreadN 40 --runMode genomeGenerate --genomeDir Genome_Dir_GRCH38 --genomeFastaFiles Genome_Index/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile Genome_Index/Homo_sapiens.GRCh38.107.gtf"__

In probedesign_processes.nf edit the params to match the location of the transcript files,genome index and specify an output folder.

To run the pipeline activate the conda environment and run: nextflow run probedesign_processes.nf 
