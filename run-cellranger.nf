#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Run CellRanger using parameters to maximize speed.
 */
process cellranger_count {
  publishDir "$output_dir", mode: 'copy'

  input:
    path transcriptome
    path fastq_path
    val count_mode
    path output_dir

  output:
    path 'cellranger-out/outs', emit: cellranger_output

  script:

    if( count_mode == 'cell' )
      """
      cellranger count --id cellranger-out --transcriptome $transcriptome --fastqs $fastq_path --nosecondary --disable-ui --nopreflight --no-bam --localcores ${task.cpus}
      """

    else if( count_mode == 'nuclear' )
      """
      cellranger count --id cellranger-out --transcriptome $transcriptome --fastqs $fastq_path --nosecondary --disable-ui --nopreflight --no-bam --localcores ${task.cpus} --include-introns
      """

    else
      error "Invalid Count Mode"



}
