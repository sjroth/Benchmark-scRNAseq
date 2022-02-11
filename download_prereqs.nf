#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Create a process that downloads the 1k PBMC test data from CellRanger.
 */
process download_testdata_1k {

  output:
    path "pbmc_1k_v3_fastqs/"
    tuple file("pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz"), file("pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz")
    tuple file("pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz"), file("pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz")
    tuple file("pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_I1_001.fastq.gz"), file("pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_I1_001.fastq.gz")

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
    path "5k_pbmc_v3_fastqs/"
    tuple file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L001_R1_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L002_R1_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L003_R1_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L004_R1_001.fastq.gz")
    tuple file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L001_R2_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L002_R2_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L003_R2_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L004_R2_001.fastq.gz")
    tuple file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L001_I1_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L002_I1_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L003_I1_001.fastq.gz"), file("5k_pbmc_v3_fastqs/5k_pbmc_v3_S1_L004_I1_001.fastq.gz")

  script:
    """
    wget https://cg.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_fastqs.tar
    tar -xvf 5k_pbmc_v3_fastqs.tar
    """
}

/*
 * Create a process to download the latest genome reference from CellRanger.
 */
process download_reference {

  output:
    path "refdata-gex-GRCh38-2020-A/"
    file "refdata-gex-GRCh38-2020-A/fasta/genome.fa"
    file "refdata-gex-GRCh38-2020-A/genes/genes.gtf"

  script:
    """
    wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
    tar -zxvf refdata-gex-GRCh38-2020-A.tar.gz
    """
}

/*
 * Create a process that will download the 10X V3 barcodes.
 */
process download_barcodes {

  output:
    file "3M-february-2018.txt"
  script:
    """
    wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz
    gunzip 3M-february-2018.txt.gz
    """
}
