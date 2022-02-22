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

    gtf_path = file.path(${gtf})
    genome_path = file.path(${genome})
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

process salmon_sel_mapping {
  input:
    path salmon_index
    path read1_files
    path read2_files
    file t2g
  output:
    path 'salmon-map', emit: salmon_map
  script:
    """
    salmon alevin -i $salmon_index -p ${task.cpus} -l IU --chromiumV3 -1 $read1_files -2 $read2_files -o salmon-map --tgMap $t2g --rad
    """
}

process salmon_sketch_mapping {
  input:
    path salmon_index
    path read1_files
    path read2_files
    file t2g
  output:
    path 'salmon-map', emit: salmon_map
  script:
    """
    salmon alevin -i $salmon_index -p ${task.cpus} -l IU --chromiumV3 -1 $read1_files -2 $read2_files -o salmon-map --tgMap $t2g --sketch
    """
}

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

process collate_rad_file_and_quant {
  input:
    path salmon_quant
    path salmon_map
    path t2g
  output:
    path 'salmon-quant-res', emit: salmon_quant_res
  script:
    """
    alevin-fry collate -t ${task.cpus} -i $salmon_quant -r $salmon_map
    alevin-fry quant -t ${task.cpus} -i $salmon_quant -o salmon-quant-res --tg-map $t2g --resolution cr-like --use-mtx
    """
}
