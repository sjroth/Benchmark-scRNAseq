#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { download_testdata_1k; download_testdata_5k; download_reference; download_barcodes } from './download_prereqs'
include { full_star_index; sparse_star_index; run_starsolo } from './run-star'

/*
 * Download prerequisites and emit their resulting downloads.
 */
workflow download_prereq_data {
  main:
    download_testdata_1k()
    download_reference()
    download_barcodes()
  emit:
    fastq_dir_1k = download_testdata_1k.out.fastq_dir
    read1_files_1k = download_testdata_1k.out.read1_files
    read2_files_1k = download_testdata_1k.out.read2_files
    index_files_1k = download_testdata_1k.out.index_files
    cellranger_reference = download_reference.out.cellranger_reference
    cellranger_genome = download_reference.out.cellranger_genome
    cellranger_gtf = download_reference.out.cellranger_gtf
    barcode_list = download_barcodes.out.barcode_list
}

/*
 * Run STAR alignment on the full genome index.
 */
 workflow star_full_index {
    take:
      genome
      gtf
      read1_files
      read2_files
      barcode_list
    main:
      full_star_index(genome, gtf)
      run_starsolo(read1_files, read2_files, barcode_list, full_star_index.out)
    emit:
      star_solo_dir = run_starsolo.out[0]
      star_log_files = run_starsolo.out[1]
 }

 /*
  * Run STAR alignment on the sparse genome index.
  */
  workflow star_sparse_index {
     take:
       genome
       gtf
       read1_files
       read2_files
       barcode_list
     main:
       sparse_star_index(genome, gtf)
       run_starsolo(read1_files, read2_files, barcode_list, sparse_star_idx.out)
     emit:
       star_solo_dir = run_starsolo.out[0]
       star_log_files = run_starsolo.out[1]
  }
