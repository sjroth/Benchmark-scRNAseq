#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { salmon_sel_mapping; salmon_sketch_mapping; salmon_quant as salmon_sel_quant; salmon_quant as salmon_sketch_quant } from './run-alevin'

workflow salmon_map_and_quant {
  take:
    salmon_index
    t2g
    read1_files
    read2_files
  main:
    salmon_sel_mapping(salmon_index,read1_files,read2_files,t2g)
    salmon_sel_quant(salmon_sel_mapping.out.salmon_map,t2g)
    salmon_sketch_mapping(salmon_index,read1_files,read2_files,t2g)
    salmon_sketch_quant(salmon_sketch_mapping.out.salmon_map,t2g)
  emit:
    salmon_sel_quant_res = salmon_sel_quant.out.salmon_quant_res
    salmon_sketch_quant_res = salmon_sketch_quant.out.salmon_quant_res
}