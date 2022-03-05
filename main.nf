#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { download_prereq_data; build_indices as build_indices_hg38;  build_indices as build_indices_mm10; get_exp as get_exp_human_1k; get_exp as get_exp_human_5k; get_exp as get_exp_human_nuclear; get_exp as get_exp_mouse } from './subworkflows'
include { kallisto_reference; run_kb_count as run_kb_count_1k; run_kb_count as run_kb_count_5k } from './run-kallisto'

outbase = "s3://fulcrumtx-users/sroth/Benchmark-scRNAseq"

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
  // Build indices.
  build_indices_hg38(download_prereq_data.out.cellranger_genome, download_prereq_data.out.cellranger_gtf)
  kallisto_reference(download_prereq_data.out.cellranger_genome, download_prereq_data.out.cellranger_gtf)

  // Get expression for cellular data.
  // 1k PBMC data.
  get_exp_human_1k(download_prereq_data.out.cellranger_reference, build_indices_hg38.out.star_idx_full, build_indices_hg38.out.star_idx_sparse, download_prereq_data.out.fastq_dir_1k, download_prereq_data.out.read1_files_1k, download_prereq_data.out.read2_files_1k, 'cell', download_prereq_data.out.barcodes_v3, 'v3', "${outbase}/1k")

  run_kb_count_1k(kallisto_reference.out.kallisto_index, kallisto_reference.out.transcripts_to_genes, download_prereq_data.out.read1_files_1k, download_prereq_data.out.read2_files_1k, "${outbase}/1k")

  // 5k PBMC data.
  get_exp_human_5k(download_prereq_data.out.cellranger_reference, build_indices_hg38.out.star_idx_full, build_indices_hg38.out.star_idx_sparse, download_prereq_data.out.fastq_dir_5k, download_prereq_data.out.read1_files_5k, download_prereq_data.out.read2_files_5k, 'cell', download_prereq_data.out.barcodes_v3, 'v3', "${outbase}/5k")

  run_kb_count_5k(kallisto_reference.out.kallisto_index, kallisto_reference.out.transcripts_to_genes, download_prereq_data.out.read1_files_5k, download_prereq_data.out.read2_files_5k, "${outbase}/5k")

  // Get expression for nuclear data.
  get_exp_human_nuclear(download_prereq_data.out.cellranger_reference, build_indices_hg38.out.star_idx_full, build_indices_hg38.out.star_idx_sparse, download_prereq_data.out.fastq_dir_nuc, download_prereq_data.out.read1_file_nuc, download_prereq_data.out.read2_file_nuc, 'nuclear', download_prereq_data.out.barcodes_v3, 'v3', "${outbase}/geo")

  // Run mouse data.
  // Build indices.
  build_indices_mm10(download_prereq_data.out.cellranger_genome_mouse, download_prereq_data.out.cellranger_gtf_mouse)

  // Get expression for mouse nuclei data.
  get_exp_mouse(download_prereq_data.out.cellranger_reference_mouse, build_indices_mm10.out.star_idx_full, build_indices_mm10.out.star_idx_sparse, download_prereq_data.out.fastq_dir_mouse_nuc, download_prereq_data.out.read1_files_mouse_nuc, download_prereq_data.out.read2_files_mouse_nuc, 'nuclear', download_prereq_data.out.barcodes_v2, 'v2', "${outbase}/mouse")
}
