#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { download_testdata_1k; download_testdata_5k; download_reference; download_barcodes } from './download_prereqs'
include { full_star_index; sparse_star_index; run_starsolo } from './run-star'
include { kallisto_reference; run_kb_count } from './run-kallisto'
include { transcriptome; transcript_to_gene; splici; remove_t2g_col; generate_salmon_index; } from './run-alevin'
include { salmon_map_and_quant; salmon_map_and_quant as salmon_quant_full_index; salmon_map_and_quant as salmon_quant_sparse_index; } from './alevin-subworkflows'

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
    run_starsolo(read1_files, read2_files, barcode_list, full_star_index.out.genome_idx)
  emit:
    star_solo_dir = run_starsolo.out.star_solo_dir
    star_log_files = run_starsolo.out.star_log_files
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
    run_starsolo(read1_files, read2_files, barcode_list, sparse_star_index.out.genome_idx)
  emit:
    star_solo_dir = run_starsolo.out.star_solo_dir
    star_log_files = run_starsolo.out.star_log_files
}

/*
 * Run kallisto pseudoalignment and quantification.
 */
workflow run_kallisto {
  take:
    genome
    gtf
    read1_files
    read2_files
  main:
    kallisto_reference(genome,gtf)
    run_kb_count(kallisto_reference.out.kallisto_index, kallisto_reference.out.transcripts_to_genes, read1_files, read2_files)
  emit:
    kallisto_output = run_kb_count.out.kallisto_output
}

/*
 * Generate transcript-to-gene mappings and salmon cDNA index. Take a genome and
 * gtf in order to generate transcripts and transcript-to-gene mapping. Use
 * transcripts to construct the salmon index.
 */
workflow salmon_cDNA_index {
  take:
    genome
    gtf
  main:
    transcriptome(genome,gtf)
    transcript_to_gene(gtf)
    generate_salmon_index(transcriptome.out.transcripts)
  emit:
    t2g = transcript_to_gene.out.t2g
    salmon_index = generate_salmon_index.out.salmon_index
}

/*
 * Generate splici transcriptome and transcript-to-gene mappings using roe.
 * Also, remove 3rd column in splici transcript-to-gene mappings for downstream
 * salmon processing.
 */
workflow splici_transcriptome {
  take:
    gtf
    genome
  main:
    splici(gtf,genome)
    remove_t2g_col(splici.out.t2g_3col)
  emit:
    transcripts = splici.out.transcripts
    t2g = remove_t2g_col.out.t2g
}

workflow salmon_cDNA {
  take:
    genome
    gtf
    read1_files
    read2_files
  main:
    salmon_cDNA_index(genome,gtf)
    salmon_map_and_quant(salmon_cDNA_index.out.salmon_index,salmon_cDNA_index.out.t2g,read1_files,read2_files)
  emit:
    salmon_sel_quant_res = salmon_map_and_quant.out.salmon_sel_quant_res
    salmon_sketch_quant_res _ salmon_map_and_quant.out.salmon_sketch_quant_res
}
