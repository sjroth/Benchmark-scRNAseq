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
  publishDir "s3://fulcrumtx-users/sroth/Benchmark-scRNAseq/", mode: "copy"

  input:
    path read1_files
    path read2_files
    file barcode_list
    path genome_idx
    val out_name

  output:
    path "${out_name}.Solo.out", emit: star_solo_dir
    tuple path("${out_name}.Log.final.out"), path("${out_name}.Log.out"), path("${out_name}.Log.progress.out"), path("${out_name}.SJ.out.tab"), emit: star_log_files

  script:
    all_r1 = "${read1_files.join(',')}"
    all_r2 = "${read2_files.join(',')}"
    """
    STAR --soloType CB_UMI_Simple --soloCBlen 16 --soloUMIlen 12 --soloCBwhitelist $barcode_list --genomeDir $genome_idx --readFilesIn $all_r2 $all_r1 --readFilesCommand zcat --outSAMtype None --runThreadN ${task.cpus} --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 --outFileNamePrefix ${out_name}.
    """
}
