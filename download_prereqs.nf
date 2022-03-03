#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Create a process that downloads the 1k PBMC test data from CellRanger.
 */
process download_testdata_1k {

  output:
    path "pbmc_1k_v3_fastqs/", emit: fastq_dir
    tuple file("pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz"), file("pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz"), emit: read1_files
    tuple file("pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz"), file("pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz"), emit: read2_files
    tuple file("pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_I1_001.fastq.gz"), file("pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_I1_001.fastq.gz"), emit: index_files

  script:
    """
    wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
    tar -xvf pbmc_1k_v3_fastqs.tar
    """
}

/*
 * Create a process that downloads the 5k PBMC test data from CellRanger.
 */
process download_testdata_5k {

  output:
    path "5k_pbmc_v3_fastqs/", emit: fastq_dir
    tuple file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L001_R1_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L002_R1_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L003_R1_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L004_R1_001.fastq.gz"), emit: read1_files
    tuple file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L001_R2_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L002_R2_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L003_R2_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L004_R2_001.fastq.gz"), emit: read2_files
    tuple file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L001_I1_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L002_I1_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L003_I1_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L004_I1_001.fastq.gz"), emit: index_files

  script:
    """
    wget https://cg.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_fastqs.tar
    tar -xvf 5k_pbmc_v3_fastqs.tar
    """
}

/*
 * Download 10X mouse E18 brain nuclei data.
 */
process download_testdata_nuclei {

  output:
    path "nuclei_900_fastqs", emit: fastq_dir
    tuple file("nuclei_900_fastqs/nuclei_900_S1_L001_R1_001.fastq.gz"), file("nuclei_900_fastqs/nuclei_900_S1_L002_R1_001.fastq.gz"), emit: read1_files
    tuple file("nuclei_900_fastqs/nuclei_900_S1_L001_R2_001.fastq.gz"), file("nuclei_900_fastqs/nuclei_900_S1_L002_R2_001.fastq.gz"), emit: read2_files
    tuple file("nuclei_900_fastqs/nuclei_900_S1_L001_I1_001.fastq.gz"), file("nuclei_900_fastqs/nuclei_900_S1_L002_I1_001.fastq.gz"), emit: index_files
  script:
    """
    wget https://cf.10xgenomics.com/samples/cell-exp/2.1.0/nuclei_900/nuclei_900_fastqs.tar
    tar -xvf nuclei_900_fastqs.tar
    """
}

/*
 * Prefetch SRA. Choose biggest file from project.
 */
process prefetch {

  output:
    path "SRR13278454/SRR13278454.sra", emit: sra_file
  script:
    """
    prefetch -p -X u -r yes -O . SRR13278454
    """
}

/*
 * Dump fastq from SRA.
 */
process fastq_dump {
  input:
    path sra_file
  output:
    tuple path("*_1.fastq"), path("*_2.fastq"), path("*_3.fastq"), emit: fastq_files
  script:
    """
    fastq-dump --split-files $sra_file
    """
}

/*
 * Gzip and rename files.
 */
process pigz {
  input:
    tuple path(fastq_1), path(fastq_2), path(fastq_3)
  output:
    path 'fastq_dir', emit: fastq_dir
    path 'fastq_dir/SRR13278454_S1_L001_R1_001.fastq.gz', emit: read1_file
    path 'fastq_dir/SRR13278454_S1_L001_R2_001.fastq.gz', emit: read2_file
    path 'fastq_dir/SRR13278454_S1_L001_I1_001.fastq.gz', emit: index_file
  script:
    """
    pigz -p ${task.cpus} $fastq_1 $fastq_2 $fastq_3
    mkdir fastq_dir
    mv ${fastq_1}.gz fastq_dir/SRR13278454_S1_L001_I1_001.fastq.gz
    mv ${fastq_2}.gz fastq_dir/SRR13278454_S1_L001_R1_001.fastq.gz
    mv ${fastq_3}.gz fastq_dir/SRR13278454_S1_L001_R2_001.fastq.gz
    """
}

/*
 * Create a process to download the latest human genome reference from CellRanger.
 */
process download_reference {

  output:
    path "refdata-gex-GRCh38-2020-A/", emit: cellranger_reference
    path "refdata-gex-GRCh38-2020-A/fasta/genome.fa", emit: cellranger_genome
    path "refdata-gex-GRCh38-2020-A/genes/genes.gtf", emit: cellranger_gtf

  script:
    """
    wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
    tar -zxvf refdata-gex-GRCh38-2020-A.tar.gz
    """
}

/*
 * Download latest mouse genome reference from CellRanger.
 */
process download_reference_mouse {

  output:
    path 'refdata-gex-mm10-2020-A/', emit: cellranger_reference
    path 'refdata-gex-mm10-2020-A/fasta/genome.fa', emit: cellranger_genome
    path 'refdata-gex-mm10-2020-A/genes/genes.gtf', emit: cellranger_gtf
  script:
    """
    wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
    tar -zxvf refdata-gex-mm10-2020-A.tar.gz
    """
}

/*
 * Create a process that will download the 10X V3 barcodes.
 */
process download_barcodes {

  output:
    path "3M-february-2018.txt", emit: barcode_list
  script:
    """
    wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz
    gunzip 3M-february-2018.txt.gz
    """
}

/*
 * Download 10X V2 barcodes.
 */
process download_barcodes_10xv2 {

  output:
    path "737K-august-2016.txt", emit: barcode_list
  script:
    """
    wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt
    """
}
