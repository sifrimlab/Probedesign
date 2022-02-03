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
params.index_files = "$baseDir/Inputs/Genome-Index"
params.kmer_length = 30
params.kmer_overlap = 20


process parse_transcripts{
    publishDir "$params.outDir", mode: 'copy'

    input:
        path transcript_file

    output:
        path "faidx_input.txt"

    script:
    """
    python3 $baseDir/src/parse_transcripts.py $transcript_file
    """

}

workflow{
  print("$baseDir")
  print(params.outDir)
  parse_transcripts("$baseDir/Inputs/transcript_file.txt")
  transcript_sequence(parse_transcripts.out)
    split_sequence(transcript_sequence.out)
    allign_probes(split_sequence.out.split_probes_for_alligning)
  // split_sequence(transcript_sequence.out.probes_for_alligning,transcript_sequence.out.complement_initial_probes)

}

process transcript_sequence{
    publishDir "$params.outDir/", mode: 'copy'

    input:
        path faidx_input
    output:
        path "probes_for_alligning.fasta", emit: probes_for_alligning
        path "complement_initial_probes.fasta", emit: complement_initial_probes

        """
        faidx  ${params.gencode_file} `cat ${faidx_input}` > probes_for_alligning.fasta
        faidx --complement ${params.gencode_file} `cat ${faidx_input}` > complement_initial_probes.fasta
        """
}

process split_sequence{
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
  publishDir "$params.outDir/", mode: 'copy'
    input:
        path split_probes_for_alligning
    output:
        path "uniquely_mapped_probes_Aligned.sam"
    script:
        """
        STAR --genomeDir ${params.index_files} --readFilesIn ${split_probes_for_alligning} \
         --outFileNamePrefix "uniquely_mapped_probes_Aligned.sam" \
         --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0
        """
}
// process select_alligned:
//     input:
//         uniquely_mapped_probes+"Aligned.out.sam"
//     output:
//         uniquely_mapped_probes_selected_sam
//         uniquely_mapped_probes_probenames
//     script:
//         '''
//         samtools view -q 10 ${input} > ${uniquely_mapped_probes_selected_sam}
//         cut -f1 ${uniquely_mapped_probes_selected_sam} > ${uniquely_mapped_probes_probenames}
//
//         samtools view -S -b -q 10 ${input} > Outputs/probes_visualizing.bam
//         samtools sort Outputs/probes_visualizing.bam  > Outputs/sorted_probes_visualizing.bam
//         samtools index Outputs/sorted_probes_visualizing.bam
//         '''
// process design_probes:
//     input:
//         split_complement_initial_probes
//         readout_file
//         barcode_file
//     output:
//         encoding_probes
//     script:
//         '''
//         python ${src_path}/probedesign.py config.yaml
//         '''
// process filter_probes:
//     input:
//         uniquely_mapped_probes_probenames
//         encoding_probes
//     output:
//         filtered_probes
//         probe_stats
//     script:
//         '''
//         python ${src_path}/filtering_probes.py config.yaml
//         '''
// process complement_sequence_probes:
//     input:
//         filtered_probes
//         uniquely_mapped_probes_probenames
//     output:
//         uniquely_mapped_probes_selected_probes_fasta
//     script:
//         '''
//         python ${src_path}/filtering_alligned_probes.py config.yaml
//         '''
// process create_blastdb:
//     input:
//         uniquely_mapped_probes_selected_probes_fasta
//     output:
//         complement_probes
//         blastdb+".nsq"
//     script:
//         '''
//         faidx --complement ${uniquely_mapped_probes_selected_probes_fasta} > ${complement_probes}
//         makeblastdb -in ${complement_probes} -dbtype nucl -parse_seqids -out ${blastdb}
//         cp ${complement_probes} ${blastdb}
//         '''
// process blast_probes:
//     input:
//         blastdb+".nsq"
//         uniquely_mapped_probes_selected_probes_fasta
//     output:
//         blast_results
//     script:
//         '''
//         blastn -db ${blastdb} -query ${uniquely_mapped_probes_selected_probes_fasta} -evalue 0.01 -out ${output}
//         '''
// process probe_results:
//     input:
//         probe_stats
//         uniquely_mapped_probes_selected_probes_fasta
//     output:
//         Final_Probe_results
//     script:
//         '''
//         python ${src_path}/results.py config.yaml
//         '''
