#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { download_prereq_data; build_indices as build_indices_hg38;  build_indices as build_indices_mm10} from './subworkflows'
include { kallisto_reference } from './run-kallisto'

/*
workflow benchmark_single_cell {
  download_prereq_data()
  star_full_index(download_prereq_data.out.cellranger_genome, download_prereq_data.out.cellranger_gtf, download_prereq_data.out.read1_files, download_prereq_data.out.read2_files, download_prereq_data.out.barcode_list)
  star_sparse_index(download_prereq_data.out.cellranger_genome, download_prereq_data.out.cellranger_gtf, download_prereq_data.out.read1_files, download_prereq_data.out.read2_files, download_prereq_data.out.barcode_list)
  cellranger_count(download_prereq_data.out.cellranger_reference, download_prereq_data.out.fastq_dir)
  run_kallisto(download_prereq_data.out.cellranger_genome, download_prereq_data.out.cellranger_gtf, download_prereq_data.out.read1_files, download_prereq_data.out.read2_files)
  salmon_cDNA(download_prereq_data.out.cellranger_genome, download_prereq_data.out.cellranger_gtf, download_prereq_data.out.read1_files, download_prereq_data.out.read2_files)
  salmon_splici(download_prereq_data.out.cellranger_gtf, download_prereq_data.out.cellranger_genome, download_prereq_data.out.read1_files, download_prereq_data.out.read2_files)
}
*/

workflow {

  // Download all prereq data.
  download_prereq_data()

  // Run human data.
  build_indices_hg38(download_prereq_data.out.cellranger_genome, download_prereq_data.out.cellranger_gtf)
  kallisto_reference(download_prereq_data.out.cellranger_genome, download_prereq_data.out.cellranger_gtf)

  // Run mouse data.
  build_indices_mm10(download_prereq_data.out.cellranger_genome_mouse, download_prereq_data.out.cellranger_gtf_mouse)
}
