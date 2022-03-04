#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { download_testdata_1k; download_testdata_5k; download_testdata_nuclei; prefetch; fastq_dump; pigz; download_reference; download_reference_mouse; download_barcodes; download_barcodes_10xv2 } from './download_prereqs'
include { full_star_index; sparse_star_index } from './run-star'
include { transcriptome; transcript_to_gene; splici; remove_t2g_col; generate_salmon_index as generate_salmon_cDNA_index; generate_salmon_index as generate_salmon_splici_index; generate_sparse_salmon_index } from './run-alevin'
include { cellranger_count } from './run-cellranger'

/*
 * Download prerequisites and emit their resulting downloads.
 */
workflow download_prereq_data {
  main:
    download_testdata_1k()
    download_testdata_5k()
    download_testdata_nuclei()

    prefetch()
    fastq_dump(prefetch.out.sra_file)
    pigz(fastq_dump.out.fastq_files)

    download_reference()
    download_reference_mouse()

    download_barcodes()
    download_barcodes_10xv2()
  emit:
    fastq_dir_1k = download_testdata_1k.out.fastq_dir
    read1_files_1k = download_testdata_1k.out.read1_files
    read2_files_1k = download_testdata_1k.out.read2_files

    fastq_dir_5k = download_testdata_5k.out.fastq_dir
    read1_files_5k = download_testdata_5k.out.read1_files
    read2_files_5k = download_testdata_5k.out.read2_files

    fastq_dir_mouse_nuc = download_testdata_nuclei.out.fastq_dir
    read1_files_mouse_nuc = download_testdata_nuclei.out.read1_files
    read2_files_mouse_nuc = download_testdata_nuclei.out.read2_files

    fastq_dir_nuc = pigz.out.fastq_dir
    read1_file_nuc = pigz.out.read1_file
    read2_file_nuc = pigz.out.read2_file

    cellranger_reference = download_reference.out.cellranger_reference
    cellranger_genome = download_reference.out.cellranger_genome
    cellranger_gtf = download_reference.out.cellranger_gtf

    cellranger_reference_mouse = download_reference_mouse.out.cellranger_reference
    cellranger_genome_mouse = download_reference_mouse.out.cellranger_genome
    cellranger_gtf_mouse = download_reference_mouse.out.cellranger_gtf

    barcodes_v3 = download_barcodes.out.barcode_list
    barcodes_v2 = download_barcodes_10xv2.out.barcode_list
}

/*
 * Build genome indices for downstream processing. Don't include kallisto cDNA reference as it is not used by both species.
 */
workflow build_indices {
  take:
    genome
    gtf
  main:
    full_star_index(genome, gtf)
    sparse_star_index(genome, gtf)

    transcriptome(genome, gtf)
    transcript_to_gene(gtf)
    generate_salmon_cDNA_index(transcriptome.out.transcripts)

    splici(gtf,genome)
    remove_t2g_col(splici.out.t2g_3col)
    generate_salmon_splici_index(splici.out.transcripts)
    generate_sparse_salmon_index(splici.out.transcripts)
  emit:
    star_idx_full = full_star_index.out.genome_idx
    star_idx_sparse = sparse_star_index.out.genome_idx

    salmon_cdna_index = generate_salmon_cDNA_index.out.salmon_index
    salmon_cdna_t2g = transcript_to_gene.out.t2g

    salmon_splici_index = generate_salmon_splici_index.out.salmon_index
    salmon_sparse_index = generate_sparse_salmon_index.out.salmon_index
    salmon_splici_t2g = remove_t2g_col.out.t2g
}

/*
 * Get expression for all workflows except kallisto cDNA because it is not used
 * for both species. Can run expression for either cellular or nuclear data.
 */
workflow get_exp {
  take:
    cellranger_reference
    fastq_path
    count_mode
    output_dir
  main:
    cellranger_count(cellranger_reference, fastq_path, count_mode, output_dir)
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
    run_starsolo(read1_files, read2_files, barcode_list, full_star_index.out.genome_idx, 'star_full_index')
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
    run_starsolo(read1_files, read2_files, barcode_list, sparse_star_index.out.genome_idx, 'star_sparse_index')
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

workflow salmon_cDNA {
  take:
    genome
    gtf
    read1_files
    read2_files
  main:
    transcriptome(genome,gtf)
    transcript_to_gene(gtf)
    generate_salmon_index(transcriptome.out.transcripts)
    salmon_map_and_quant(generate_salmon_index.out.salmon_index, transcript_to_gene.out.t2g, read1_files, read2_files, 'salmon-cDNA')
  emit:
    salmon_sel_quant_res = salmon_map_and_quant.out.salmon_sel_quant_res
    salmon_sketch_quant_res = salmon_map_and_quant.out.salmon_sketch_quant_res
}

workflow salmon_splici {
  take:
    gtf
    genome
    read1_files
    read2_files
  main:
    splici(gtf,genome)
    remove_t2g_col(splici.out.t2g_3col)
    generate_salmon_index(splici.out.transcripts)
    generate_sparse_salmon_index(splici.out.transcripts)
    salmon_quant_full_index(generate_salmon_index.out.salmon_index, remove_t2g_col.out.t2g, read1_files, read2_files, 'salmon-splici-full')
    salmon_quant_sparse_index(generate_sparse_salmon_index.out.salmon_index, remove_t2g_col.out.t2g, read1_files, read2_files, 'salmon-splici-sparse')
  emit:
    salmon_sel_full_quant_res = salmon_quant_full_index.out.salmon_sel_quant_res
    salmon_sketch_full_quant_res = salmon_quant_full_index.out.salmon_sketch_quant_res
    salmon_sel_sparse_quant_res = salmon_quant_sparse_index.out.salmon_sel_quant_res
    salmon_sketch_sparse_quant_res = salmon_quant_sparse_index.out.salmon_sketch_quant_res
}
