#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Use gffread to generate a transcripts file from the CellRanger genome and gtf.
 */
process transcriptome {
  input:
    file genome
    file gtf
  output:
    path 'transcripts.fa', emit: transcripts
  script:
    """
    gffread -w transcripts.fa -g $genome $gtf
    """
}

/*
 * Generate splici genome and transcript-to-gene mapping for using in splici
 * index and downstream processing.
 */
process splici {
  input:
    file gtf
    file genome
  output:
    path 'transcriptome_splici_fl86.fa', emit: transcripts
    path 'transcriptome_splici_fl86_t2g_3col.tsv', emit: t2g_3col
  script:
    """
    #!/usr/bin/env Rscript

    library("roe")

    gtf_path = file.path('${gtf}')
    genome_path = file.path('${genome}')
    read_length = 91
    flank_trim_length = 5

    make_splici_txome(gtf_path=gtf_path,
                      genome_path=genome_path,
                      read_length=read_length,
                      flank_trim_length=flank_trim_length,
                      output_dir=".",
                      file_name="transcriptome_splici")
    """
}

/*
 * Remove third column of splici transcript-to-gene mapping.
 */
process remove_t2g_col {
  input:
    file t2g_3col
  output:
    path 't2g.tsv', emit: t2g
  script:
    """
    cut -f1,2 $t2g_3col > t2g.tsv
    """
}

/*
 * Generate transcript to gene mapping.
 */
process transcript_to_gene {
  input:
    file gtf
  output:
    path 't2g.tsv', emit: t2g
  shell:
    template('t2g.txt')
}

/*
 * Generate a salmon index with default parameters.
 */
process generate_salmon_index {
  input:
    file transcripts
  output:
    path 'salmon-index', emit: salmon_index
  script:
    """
    salmon index -t $transcripts -i salmon-index -p ${task.cpus}
    """
}

/*
 * Generate a sparse salmon index.
 */
process generate_sparse_salmon_index {
  input:
    file transcripts
  output:
    path 'salmon-index', emit: salmon_index
  script:
    """
    salmon index -t $transcripts -i salmon-index -p ${task.cpus} --sparse
    """
}

/*
 * Perform selective mapping.
 */
process salmon_sel_mapping {
  input:
    path salmon_index
    path read1_files
    path read2_files
    file t2g
    val chemistry
  output:
    path 'salmon-map', emit: salmon_map
  script:

    if( chemistry == 'v3' )
      chem_cmd = '--chromiumV3'
    else if( chemistry == 'v2' )
      chem_cmd = '--chromium'
    else
      error "Invalid 10X Chemistry"

    """
    salmon alevin -i $salmon_index -p ${task.cpus} -l IU $chem_cmd -1 $read1_files -2 $read2_files -o salmon-map --tgMap $t2g --rad
    """
}

/*
 * Perform sketch mapping.
 */
process salmon_sketch_mapping {
  input:
    path salmon_index
    path read1_files
    path read2_files
    file t2g
    val chemistry
  output:
    path 'salmon-map', emit: salmon_map
  script:

    if( chemistry == 'v3' )
      chem_cmd = '--chromiumV3'
    else if( chemistry == 'v2' )
      chem_cmd = '--chromium'
    else
      error "Invalid 10X Chemistry"

    """
    salmon alevin -i $salmon_index -p ${task.cpus} -l IU $chem_cmd -1 $read1_files -2 $read2_files -o salmon-map --tgMap $t2g --sketch
    """
}

/*
 * Get list of called cell barcodes.
 */
process generate_permit_list {
  input:
    path salmon_map
  output:
    path 'salmon-quant', emit: salmon_quant
  script:
    """
    alevin-fry generate-permit-list -d fw -k -i $salmon_map -o salmon-quant
    """
}

/*
 * Perform quantification.
 */
process collate_rad_file_and_quant {

  input:
    path salmon_quant
    path salmon_map
    file t2g
  output:
    path "alevin-out/alevin/quants_mat.mtx", emit: alevin_mtx
    path "alevin-out/alevin/quants_mat_cols.txt", emit: alevin_features
    path "alevin-out/alevin/quants_mat_rows.txt", emit: alevin_barcodes
  script:
    """
    alevin-fry collate -t ${task.cpus} -i $salmon_quant -r $salmon_map
    alevin-fry quant -t ${task.cpus} -i $salmon_quant -o alevin-out --tg-map $t2g --resolution cr-like --use-mtx
    """
}

process format_alevin_output {
  publishDir "s3://fulcrumtx-users/sroth/Benchmark-scRNAseq/", mode: "copy"
  
  input:
    path alevin_mtx
    path alevin_features
    path alevin_barcodes
    val output
  output:
    path "alevin-${output}", emit: alevin_outdir
    path "alevin-${output}/matrix.mtx.gz", emit: alevin_out_mtx
    path "alevin-${output}/features.tsv.gz", emit: alevin_out_features
    path "alevin-${output}/barcodes.tsv.gz", emit: alevin_out_barcodes
  script:
    """
    mkdir alevin-${output}
    cp $alevin_mtx alevin-${output}/matrix.mtx
    cp $alevin_features alevin-${output}/features.tsv
    cp $alevin_barcodes alevin-${output}/barcodes.tsv
    gzip alevin-${output}
    """
}

/*
 * Perform salmom quantification after salmon mapping.
 */
workflow salmon_quant {
  take:
    salmon_map
    t2g
    output
  main:
    generate_permit_list(salmon_map)
    collate_rad_file_and_quant(generate_permit_list.out.salmon_quant, salmon_map, t2g)
    format_alevin_output(collate_rad_file_and_quant.out.alevin_mtx, collate_rad_file_and_quant.out.alevin_features, collate_rad_file_and_quant.out.alevin_barcodes, output)
  emit:
    alevin_outdir = format_alevin_output.out.alevin_outdir
}
