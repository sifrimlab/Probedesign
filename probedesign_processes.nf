nextflow.enable.dsl=2

// process all:
//     input:
//         probes_for_alligning
//         complement_initial_probes
//         split_probes_for_alligning
//         split_complement_initial_probes
//         uniquely_mapped_probes+"Aligned.out.sam"
//         uniquely_mapped_probes_probenames
//         uniquely_mapped_probes_selected_sam
//         encoding_probes
//         filtered_probes
//         probe_stats
//         uniquely_mapped_probes_selected_probes_fasta
//         complement_probes
//         blastdb+".nsq"
//         blast_results
//         Final_Probe_results

params.outDir = "$baseDir/nf_output"
params.gencode_file = "$baseDir/Inputs/gencode.v37.transcripts.fa"
// params.transcript_file = "$baseDir/Inputs/transcript_file_2.txt"
params.transcript_file = "$baseDir/Inputs/transcript_file.txt"
params.index_files = "/media/yob/Genome_Dir_GRCH38"
params.readout_file = "Inputs/Readout_probes_information.csv"
params.barcode_file = "Inputs/barcodes_merfish.csv"
params.kmer_length = 30
params.kmer_overlap = 20


process parse_transcripts{
    cpus 8
    publishDir "$params.outDir", mode: 'copy'

    input:
        path transcript_file

    output:
        path "probes_for_alligning.fasta", emit: probes_for_alligning
        path "complement_initial_probes.fasta", emit: complement_initial_probes

    script:
    """
    python3 $baseDir/src/parse_transcripts.py $params.transcript_file $params.gencode_file
    """
}

process split_sequence{
    cpus 8
  publishDir "$params.outDir/", mode: 'copy'
    input:
        path probes_for_alligning
        path complement_initial_probes
    output:
        path "split_probes_for_alligning.fasta", emit: split_probes_for_alligning
        path "split_complement_initial_probes.fasta", emit: split_complement_initial_probes

    script:
        """
        pyfasta split -n1 -k${params.kmer_length} -o${params.kmer_overlap} ${probes_for_alligning}
        mv probes_for_alligning.split.${params.kmer_length}mer.${params.kmer_overlap}overlap.fasta split_probes_for_alligning.fasta
        pyfasta split -n1 -k${params.kmer_length} -o${params.kmer_overlap} ${complement_initial_probes}
        mv complement_initial_probes.split.${params.kmer_length}mer.${params.kmer_overlap}overlap.fasta split_complement_initial_probes.fasta
        """
}

process allign_probes{
    cpus 8
  publishDir "$params.outDir/", mode: 'copy'
    input:
        path split_probes_for_alligning
    output:
        path "uniquely_mapped_probes_Aligned.out.sam", emit: uniquely_mapped_probes_Aligned
    script:
        """
        STAR --genomeDir ${params.index_files} --readFilesIn ${split_probes_for_alligning} \
         --outFileNamePrefix "uniquely_mapped_probes_" \
         --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0
        """
}

process select_alligned{
    cpus 8
    publishDir "$params.outDir/", mode: 'copy'
    input:
      path uniquely_mapped_probes_Aligned
    output:
      path  "uniquely_mapped_probes_selected.sam"
      path  "uniquely_mapped_probes_probenames.txt", emit: uniquely_mapped_probes_probenames
      path  "sorted_probes_visualizing.bam"
    script:
        """
        samtools view -q 10 ${uniquely_mapped_probes_Aligned} > uniquely_mapped_probes_selected.sam
        cut -f1 uniquely_mapped_probes_selected.sam > uniquely_mapped_probes_probenames.txt

        samtools view -S -b -q 10 ${uniquely_mapped_probes_Aligned} > probes_visualizing.bam
        samtools sort probes_visualizing.bam  > sorted_probes_visualizing.bam
        samtools index sorted_probes_visualizing.bam
        """

}

process design_probes{
    cpus 8
    publishDir "$params.outDir/", mode: 'copy'
    input:
        path split_complement_initial_probes
    output:
        path "encoding_probes.fasta", emit: encoding_probes
    script:
        """
        python3 $baseDir/src/probedesign.py $baseDir/config.yaml $baseDir $baseDir/nf_output $split_complement_initial_probes
        """
}

process filter_probes{
    cpus 8
    publishDir "$params.outDir/", mode: 'copy'
    input:
        path encoding_probes
        path uniquely_mapped_probes_probenames
    output:
        path "filtered_probes.fasta", emit: filtered_probes
        path "probe_stats.csv", emit: probe_stats
    script:
        """
        python $baseDir/src/filtering_probes.py $baseDir/config.yaml $baseDir $baseDir/nf_output  $uniquely_mapped_probes_probenames
        """
}

process complement_sequence_probes{
    cpus 8
    publishDir "$params.outDir/", mode: 'copy'
    input:
        path filtered_probes
        path uniquely_mapped_probes_probenames
    output:
        path "uniquely_mapped_probes_filtered_final.fasta", emit: uniquely_mapped_probes_filtered_final
    script:
        """
        python $baseDir/src/filtering_alligned_probes.py $baseDir/config.yaml $baseDir $baseDir/nf_output
        """
}
process create_blastdb{
    cpus 8
    publishDir "$params.outDir/", mode: 'copy'
    input:
        path uniquely_mapped_probes_filtered_final
    output:
         path "complement_probes*"
//        path "complement_probes.fasta"
//        path "complement_probes_blastdb.nsq", emit:complement_probes_blastdb
    script:
        """
        mkdir -p $baseDir/blastdb
        faidx --complement ${uniquely_mapped_probes_filtered_final} > complement_probes.fasta
        makeblastdb -in complement_probes.fasta -dbtype nucl -parse_seqids -out complement_probes_blastdb
        """
}

process blast_probes{
    cpus 8
    publishDir "$params.outDir/", mode: 'copy'
    input:
       path complement_probes_blastdb
       path uniquely_mapped_probes_filtered_final
    output:
       path "blast_results.out"
    script:
        """
        blastn -db complement_probes_blastdb -query ${uniquely_mapped_probes_filtered_final} -evalue 0.01 -out blast_results.out
        """
}

process probe_results{
    cpus 8
    publishDir "$params.outDir/", mode: 'copy'
    input:
      path probe_stats
      path uniquely_mapped_probes_filtered_final
    output:
      path "Final_Probe_results.csv"
    script:
      """
      python $baseDir/src/results.py $baseDir/config.yaml $baseDir $baseDir/nf_output
      """
}

workflow{
  parse_transcripts(params.transcript_file)
  split_sequence(parse_transcripts.out)
  allign_probes(split_sequence.out.split_probes_for_alligning)
  select_alligned(allign_probes.out)
  design_probes(split_sequence.out.split_complement_initial_probes)
  filter_probes(design_probes.out,select_alligned.out.uniquely_mapped_probes_probenames)
  complement_sequence_probes(filter_probes.out.filtered_probes,select_alligned.out.uniquely_mapped_probes_probenames)
  create_blastdb(complement_sequence_probes.out)
  blast_probes(create_blastdb.out.collect(),complement_sequence_probes.out)
  probe_results(filter_probes.out.probe_stats,complement_sequence_probes.out)
}
