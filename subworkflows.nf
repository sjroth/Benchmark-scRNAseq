#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { download_testdata_1k, download_testdata_5k, download_reference, download_barcodes } from './download_prereqs.nf'
include { full_star_index, sparse_star_index, run_starsolo } from './run-star.nf'

/*
 * Download prerequisites and emit their resulting downloads.
 */
workflow download_prereqs {
  main:
    download_testdata_1k()
    download_reference()
    download_barcodes()
  emit:
    1k_fastq_dir = download_testdata_1k.out[0]
    1k_read1_files = download_testdata_1k.out[1]
    1k_read2_files = download_testdata_1k.out[2]
    1k_index_files = download_testdata_1k.out[3]
    cellranger_reference = download_reference.out[0]
    cellranger_genome = download_reference.out[1]
    cellranger_gtf = download_reference.out[2]
    barcode_list = download_barcodes.out
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
       sparse_star_idx(genome, gtf)
       run_starsolo(read1_files, read2_files, barcode_list, sparse_star_idx.out)
     emit:
       star_solo_dir = run_starsolo.out[0]
       star_log_files = run_starsolo.out[1]
  }
