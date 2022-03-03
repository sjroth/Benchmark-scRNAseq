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
      count_cmd = ''

    else if( mode == 'nuclear' )
      count_cmd = '--include-introns'

    else
      error "Invalid Count Mode"

    """
    cellranger count --id cellranger-out --transcriptome $transcriptome --fastqs $fastq_path --nosecondary --disable-ui --nopreflight --no-bam --localcores ${task.cpus} $count_cmd
    """

}
