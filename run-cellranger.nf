#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Run CellRanger using parameters to maximize speed.
 */
process cellranger_count {
  publishDir "s3://fulcrumtx-users/sroth/Benchmark-scRNAseq/", mode: 'copy'

  input:
    path transcriptome
    path fastq_path
    val count_mode
    val output

  output:
    path "${output}/cellranger", emit: cellranger_output

  script:

    if( count_mode == 'cell' )
      count_cmd = ''

    else if( count_mode == 'nuclear' )
      count_cmd = '--include-introns'

    else
      error "Invalid Count Mode"

    """
    cellranger count --id out --transcriptome $transcriptome --fastqs $fastq_path --nosecondary --disable-ui --nopreflight --no-bam --localcores ${task.cpus} $count_cmd
    mkdir ${output}
    mv out/outs ${output}/cellranger
    """

}
