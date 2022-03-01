#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { download_prereq_data; star_full_index; star_sparse_index; run_kallisto; salmon_cDNA; salmon_splici } from './subworkflows'
include { cellranger_count } from './run-cellranger'

workflow {
  download_prereq_data()
  star_full_index(download_prereq_data.out.cellranger_genome, download_prereq_data.out.cellranger_gtf, download_prereq_data.out.read1_files, download_prereq_data.out.read2_files, download_prereq_data.out.barcode_list)
  star_sparse_index(download_prereq_data.out.cellranger_genome, download_prereq_data.out.cellranger_gtf, download_prereq_data.out.read1_files, download_prereq_data.out.read2_files, download_prereq_data.out.barcode_list)
  cellranger_count(download_prereq_data.out.cellranger_reference, download_prereq_data.out.fastq_dir)
  run_kallisto(download_prereq_data.out.cellranger_genome, download_prereq_data.out.cellranger_gtf, download_prereq_data.out.read1_files, download_prereq_data.out.read2_files)
  salmon_cDNA(download_prereq_data.out.cellranger_genome, download_prereq_data.out.cellranger_gtf, download_prereq_data.out.read1_files, download_prereq_data.out.read2_files)
  salmon_splici(download_prereq_data.out.cellranger_gtf, download_prereq_data.out.cellranger_genome, download_prereq_data.out.read1_files, download_prereq_data.out.read2_files)
}
/*

workflow {
  download_reference()
  kallisto_reference_nuclear(download_reference.out.cellranger_genome,download_reference.out.cellranger_gtf)
}
*/
