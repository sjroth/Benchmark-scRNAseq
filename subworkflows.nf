#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { download_testdata_1k; download_testdata_5k; download_testdata_nuclei; prefetch; fastq_dump; pigz; download_reference; download_reference_mouse; download_barcodes; download_barcodes_10xv2 } from './download_prereqs'
include { full_star_index; sparse_star_index; run_starsolo as run_starsolo_full; run_starsolo as run_starsolo_sparse; format_star_output as format_star_output_full; format_star_output as format_star_output_sparse } from './run-star'
include { transcriptome; transcript_to_gene; splici; remove_t2g_col; generate_salmon_index as generate_salmon_cDNA_index; generate_salmon_index as generate_salmon_splici_index; generate_sparse_salmon_index } from './run-alevin'
include { cellranger_count } from './run-cellranger'
include { salmon_map_and_quant as alevin_quant_cDNA; salmon_map_and_quant as alevin_quant_splici_full; salmon_map_and_quant as alevin_quant_splici_sparse } from './alevin-subworkflows'

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
    // Genome indices.
    cellranger_reference
    star_idx_full
    star_idx_sparse
    salmon_cdna_index
    salmon_splici_index
    salmon_sparse_index

    // Data paths
    fastq_path
    read1_files
    read2_files

    // 10X options
    count_mode
    barcode_list
    chemistry

    // Alevin transcript-to-gene mappings
    salmon_cdna_t2g
    salmon_splici_t2g

    // Output
    output

  main:

    cellranger_count(cellranger_reference, fastq_path, count_mode, output)

    run_starsolo_full(read1_files, read2_files, barcode_list, star_idx_full, count_mode, chemistry)
    format_star_output_full(run_starsolo_full.out.star_solo_mtx, run_starsolo_full.out.star_solo_features, run_starsolo_full.out.star_solo_barcodes, "${output}/STARsolo-full")

    run_starsolo_sparse(read1_files, read2_files, barcode_list, star_idx_sparse, count_mode, chemistry)
    format_star_output_sparse(run_starsolo_sparse.out.star_solo_mtx, run_starsolo_sparse.out.star_solo_features, run_starsolo_sparse.out.star_solo_barcodes, "${output}/STARsolo-sparse")

    alevin_quant_cDNA(salmon_cdna_index, salmon_cdna_t2g, read1_files, read2_files, chemistry, "${output}/alevin-cDNA")

    alevin_quant_splici_full(salmon_splici_index, salmon_splici_t2g, read1_files, read2_files, chemistry, "${output}/alevin-splici-full")
    alevin_quant_splici_sparse(salmon_sparse_index, salmon_splici_t2g, read1_files, read2_files, chemistry, "${output}/alevin-splici-sparse")

  emit:
    cellranger_out = cellranger_count.out.cellranger_output

    star_solo_outdir_full = format_star_output_full.out.star_solo_outdir
    star_solo_outdir_sparse = format_star_output_sparse.out.star_solo_outdir

    alevin_cDNA_sel_outdir = alevin_quant_cDNA.out.salmon_sel_outdir
    alevin_cDNA_sketch_quant = alevin_quant_cDNA.out.salmon_sketch_outdir

    alevin_splici_full_sel_quant = alevin_quant_splici_full.out.salmon_sel_outdir
    alevin_splici_full_sketch_quant = alevin_quant_splici_full.out.salmon_sketch_outdir

    alevin_splici_sparse_sel_quant = alevin_quant_splici_sparse.out.salmon_sel_outdir
    alevin_splici_sparse_sketch_quant = alevin_quant_splici_sparse.out.salmon_sketch_outdir
}
