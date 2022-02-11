#! /usr/bin/env nextflow

include { run_starsolo } from './run-star.nf'
include { download_testdata_1k } from './download_prereqs.nf'

nextflow.enable.dsl = 2

workflow {
  download_testdata_1k()
  run_starsolo( download_testdata_1k.out[1], download_testdata_1k.out[2] )
}
