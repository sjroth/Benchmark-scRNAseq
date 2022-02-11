#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { download_prereq_data } from './subworkflows'

workflow {
  download_prereq_data()
}
