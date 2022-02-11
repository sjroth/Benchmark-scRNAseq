#! /usr/bin/env nextflow

include { download_prereqs, star_full_index, star_sparse_index } from './subworkflows.nf'

nextflow.enable.dsl = 2

workflow {
  download_prereqs()
  star_full_index(download_prereqs.out.cellranger_genome, download_prereqs.out.cellranger_gtf, download_prereqs.out.1k_read1_files, download_prereqs.out.1k_read2_files, download_prereqs.out.barcode_list)
  star_sparse_index(download_prereqs.out.cellranger_genome, download_prereqs.out.cellranger_gtf, download_prereqs.out.1k_read1_files, download_prereqs.out.1k_read2_files, download_prereqs.out.barcode_list)
}
