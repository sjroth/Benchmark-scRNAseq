#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Construct a Kallisto reference.
 */
process kallisto_reference {
  input:
    file genome
    file gtf
  output:
    path 'transcriptome.idx', emit: kallisto_index
    path 'transcripts_to_genes.txt', emit: transcripts_to_genes
    path 'cdna.fa', emit: cdna
  script:
    """
    kb ref -i transcriptome.idx -g transcripts_to_genes.txt -f1 cdna.fa $genome $gtf
    """
}

/*
 * Quantify using Kallisto.
 */
process run_kb_count {
  input:
    path kallisto_index
    path transcripts_to_genes
    path read1_files
    path read2_files
  output:
    path 'kallisto-test', emit: kallisto_output
  script:
    read1_lst = read1_files.toList()
    read2_lst = read2_files.toList()
    all_fastq = [read1_lst, read2_lst].transpose().flatten()

    """
    kb count -i $kallisto_index -g $transcripts_to_genes -x 10XV3 -o kallisto-test -t ${task.cpus} ${all_fastq.join(' ')}
    """

}
