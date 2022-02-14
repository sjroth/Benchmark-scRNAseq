#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Run CellRanger using parameters to maximize speed.
 */
process cellranger_count {

  input:
    path transcriptome
    path fastq_path

  output:
    path 'cellranger-test/outs', emit: cellranger_output

  script:
    """
    cellranger count --id cellranger-test --transcriptome $transcriptome --fastqs $fastq_path --nosecondary --disable-ui --nopreflight --no-bam --localcores ${task.cpus}
    """

}
