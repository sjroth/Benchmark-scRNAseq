#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { download_prereq_data, star_full_index } from './subworkflows'

workflow {
  download_prereq_data()
  star_full_index(download_prereq_data.out.cellranger_genome, download_prereq_data.out.cellranger_gtf,download_prereq_data.out.read1_files_1k,download_prereq_data.out.read2_files_1k,download_prereq_data.out.barcode_list)
}
