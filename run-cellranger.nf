#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Run CellRanger using parameters to maximize speed.
 */
process cellranger_count {
  publishDir 's3://fulcrumtx-users/sroth/Benchmark-scRNAseq/5k/', mode: 'copy'

  input:
    path transcriptome
    path fastq_path

  output:
    path 'cellranger-out/outs', emit: cellranger_output

  script:
    """
    cellranger count --id cellranger-out --transcriptome $transcriptome --fastqs $fastq_path --nosecondary --disable-ui --nopreflight --no-bam --localcores ${task.cpus}
    """

}
