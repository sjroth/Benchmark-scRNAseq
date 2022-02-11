#! /usr/bin/env nextflow

include { download_prereqs, star_full_index, star_sparse_index } from './subworkflows.nf'

nextflow.enable.dsl = 2

workflow {
  download_prereqs()
  
}
