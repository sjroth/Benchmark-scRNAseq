#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Construct a Kallisto reference for standard scRNA-seq.
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
    path "kallisto-output/counts_unfiltered/cells_x_genes.mtx", emit: kallisto_mtx
    path "kallisto-output/counts_unfiltered/cells_x_genes.genes.txt", emit: kallisto_features
    path "kallisto-output/counts_unfiltered/cells_x_genes.barcodes.txt", emit: kallisto_barcodes
  script:
    read1_lst = read1_files.toList()
    read2_lst = read2_files.toList()
    all_fastq = [read1_lst, read2_lst].transpose().flatten()

    """
    kb count -i $kallisto_index -g $transcripts_to_genes -x 10XV3 -o kallisto-output -t ${task.cpus} ${all_fastq.join(' ')}
    """

}

process format_kb_output {
  publishDir "s3://fulcrumtx-users/sroth/Benchmark-scRNAseq/", mode: 'copy'

  input:
    path kallisto_mtx
    path kallisto_features
    path kallisto_barcodes
    val output
  output:
    path "kallisto-${output}", emit: kallisto_outdir
    path "kallisto-${output}/matrix.mtx.gz", emit: kallisto_out_mtx
    path "kallisto-${output}/features.tsv.gz", emit: kallisto_out_features
    path "kallisto-${output}/barcodes.tsv.gz", emit: kallisto_out_barcodes
  script:
    """
    mkdir kallisto-${output}
    cp $kallisto_mtx kallisto-${output}/matrix.mtx
    cp $kallisto_features kallisto-${output}/features.tsv
    cp $kallisto_barcodes kallisto-${output}/barcodes.tsv
    gzip kallisto-${output}/*
    """
}
