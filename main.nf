#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/circrnasponging
========================================================================================
 nf-core/circrnasponging Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/circrnasponging
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/circrnasponging --input '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --input [file]                  Path to input data (must be surrounded with quotes)
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, docker, singularity, test, awsbatch, <institute> and more
    Options:
      --genome [str]                  Name of iGenomes reference
      --single_end [bool]             Specifies that the input is single-end reads

    References                        If not specified in the configuration file or you wish to overwrite any of the references
      --fasta [file]                  Path to fasta reference

    Other options:
      --outdir [file]                 The output directory where the results will be saved
    """.stripIndent() 
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * Create a channel for input read files
 */

Channel
        .fromFilePairs(params.input, size: params.single_end ? 1 : 2)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
        .set{ch_totalRNA_reads}

process STAR {
    tag "$name"
    //label 'process_medium'
    publishDir "${params.out_dir}/samples/${sampleID}/circRNA_detection/", mode: params.publish_dir_mode
   //    saveAs: { filename ->
   //                   filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
   //             }

    input:
    set val(sampleID), file(reads) from ch_totalRNA_reads

    output:
    tuple val(sampleID), file("Chimeric.out.junction") into chimeric_junction_files

    script:
    """
    STAR --chimSegmentMin 10 --runThreadN 10 --genomeDir $params.STAR_index --readFilesCommand zcat --readFilesIn $reads
    """
}

