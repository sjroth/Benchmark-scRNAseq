#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { download_prereq_data } from './subworkflows'
include { download_testdata_1k } from './download_prereqs'

workflow {
  download_prereq_data()
}
