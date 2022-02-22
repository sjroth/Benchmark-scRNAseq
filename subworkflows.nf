#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { download_testdata_1k; download_testdata_5k; download_reference; download_barcodes } from './download_prereqs'
include { full_star_index; sparse_star_index; run_starsolo } from './run-star'
include { kallisto_reference; run_kb_count } from './run-kallisto'
include { transcriptome; transcript_to_gene; generate_salmon_index; salmon_sel_mapping; salmon_sketch_mapping ; generate_permit_list; collate_rad_file_and_quant } from './run-alevin'

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

workflow salmon_cDNA_sel {
  take:
    salmon_index
    read1_files
    read2_files
    t2g
  main:
    salmon_sel_mapping(salmon_index,read1_files,read2_files,t2g)
    generate_permit_list(salmon_sel_mapping.out.salmon_map)
    collate_rad_file_and_quant(generate_permit_list.out.salmon_quant,salmon_sel_mapping.out.salmon_map,t2g)
  emit:
    salmon_out = collate_rad_file_and_quant.out.salmon_quant_res
}

workflow salmon_cDNA_sketch {
  take:
    salmon_index
    read1_files
    read2_files
    t2g
  main:
    salmon_sketch_mapping(salmon_index,read1_files,read2_files,t2g)
    generate_permit_list(salmon_sketch_mapping.out.salmon_map)
    collate_rad_file_and_quant(generate_permit_list.out.salmon_quant,salmon_sketch_mapping.out.salmon_map,t2g)
  emit:
    salmon_out = collate_rad_file_and_quant.out.salmon_quant_res
}

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
