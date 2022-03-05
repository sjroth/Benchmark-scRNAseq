#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Create a process that constructs a full STAR index.
 */
process full_star_index {

  input:
    file genome
    file gtf

  output:
    path 'genome_idx', emit: genome_idx

  script:
    """
    mkdir genome_idx
    STAR --runMode genomeGenerate --runThreadN ${task.cpus} --genomeDir genome_idx --genomeFastaFiles $genome --sjdbGTFfile $gtf
    """

}

/*
 * Create a process that constructs a sparse STAR index to best mimic CellRanger.
 */
process sparse_star_index {

  input:
    file genome
    file gtf

  output:
    path 'genome_idx', emit: genome_idx

  script:
    """
    mkdir genome_idx
    STAR --runMode genomeGenerate --runThreadN ${task.cpus} --genomeDir genome_idx --genomeFastaFiles $genome --sjdbGTFfile $gtf --genomeSAsparseD 3
    """

}

/*
 * Create a process that performs STARsolo. Use settings to best mimic CellRanger.
 */
process run_starsolo {
  publishDir "$output_dir", mode: "copy"

  input:
    path read1_files
    path read2_files
    file barcode_list
    path genome_idx
    val count_mode
    val chemistry
    path output_dir
    val out_name

  output:
    path "${out_name}.Solo.out", emit: star_solo_dir
    tuple path("${out_name}.Log.final.out"), path("${out_name}.Log.out"), path("${out_name}.Log.progress.out"), path("${out_name}.SJ.out.tab"), emit: star_log_files

  script:
    all_r1 = "${read1_files.join(',')}"
    all_r2 = "${read2_files.join(',')}"

    if( count_mode == 'cell' )
      count_cmd = ''
    else if( count_mode == 'nuclear' )
      count_cmd = '--soloFeatures GeneFull'
    else
      error "Invalid Count Mode"

    if( chemistry == 'v3' )
      umi_len = '12'
    else if( chemistry == 'v2' )
      umi_len = '10'
    else
      error "Invalid 10X Chemistry"

    """
    STAR --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen $umi_len --soloBarcodeReadLength 0 --soloCBwhitelist $barcode_list --genomeDir $genome_idx --readFilesIn $all_r2 $all_r1 $count_cmd --readFilesCommand zcat --outSAMtype None --runThreadN ${task.cpus} --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 --outFileNamePrefix ${out_name}.
    """
}
