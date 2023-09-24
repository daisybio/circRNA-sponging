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

nextflow.enable.dsl=2

def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/circrnasponging 

    Mandatory arguments:
      --samplesheet [file]		Path to samplesheet (must be surrounded with quotes)
      --outdir [file]			The output directory where the results will be saved
      --species [str]			Species name in 3 letter code (hsa for human, mmu for mouse)
      --genome [str]    Genome version that will be used for mapping e.g. for human hg38
      --miRNA_adapter [str] 		miRNA adapter used for trimming
      -profile [str]           	 	Configuration profile to use. Can use multiple (comma separated)
                                      	Available: conda, docker, singularity, test, awsbatch, <institute> and more

      --fasta [file] 			Path to genome fasta (must be surrounded with quotes)
      --gtf [file]			Path to gtf file (must be surrounded with quotes)
      --bed12 [file]		Path to gene annotation (must be surrounded with quotes)
      --miRNA_fasta [file]		Path to mature miRNA fasta (must be surrounded with quotes)
      --miRNA_related_fasta [file]	Path to mature miRNA fasta of related species (must be surrounded with quotes)
      --hairpin_fasta [file]		Path to miRNA hairpin fasta (must be surrounded with quotes)
      --transcriptome [file]        Path to transcriptome fasta (must be surrounded with quotes)

    Options:
      --miRNA_raw_counts [file]		Path to tabulated raw miRNA counts (must be surrounded with quotes)
      --single_end [bool]            	Specifies that the total RNA input is single-end reads
      --sample_group	[file]		File specifying partitioning of samples into groups (must be surrounded with quotes)
      --read_threshold [real]		Positive. Read counts under this threshold are considered to be low expressed
      --sample_percentage [real]	Between 0 and 1. Minimum percentage of samples that should have no low expression
      --circRNA_only [bool]  		Run only circRNA analysis, don't run miRNA analysis
      --splitter [int]          Number of fasta circRNA entries to use for each miRNA target prediction run
      --quantification [bool]       Use psirc-quant quantified circRNA expression for downstream analysis (default: true)
      --database_annotation [bool]  Annotate circRNA hits with circBase data
        --annotated_only [bool]      Only use circRNAs that could be annotated by CircBase for downstream analysis (default: false)
        --offline_circ_db [file]      File containing downloaded circBase entries for offline access to the database
      --differential_expression [bool]  Enable differential expression analysis using DESeq2 on all given RNA-seq data and circRNA only
      miRNA binding site predictions:
        --tarpmir [bool]   Use TarPmiR
            --model [pkl]       Path to a specific tarpmir model as pickle (default model is Human_sklearn_0.22.pkl, located in data directory)
            --p     [real]      Double between 0 and 1 specifing cutoff for tarpmir (default: 0.8)
            --threads [int]     Number of threads to use for tarpmir (default: 2)
        --pita [bool]  Use PITA
            --pita_l [str]      Search for seed lengths of num1,...,num2 to the MicroRNA (default: "6-8")
            --pita_gu [str]     Lengths for which G:U wobbles are allowed and number of allowed wobbles.
                                    Format of nums: <length;num G:U>,<length;num G:U>,... (default: "6;0,7;1,8;1")
            --pita_m [str]      Lengths for which mismatches are allowed and number of allowed mismatches
                                    Format of nums: <length;num mismatches>,<length;num mismatches>,...
                                    (default: "6;0,7;0,8;1")
      --sponge [bool]   Wheather to perform SPONGE ceRNA network analysis    
        --target_scan_symbols [file]    Path to target scan symbols matrix in tsv and SPONGE format (rows=ENSG,cols=miRNA,data=counts)
        --tpm [bool]     Use log(TPM + 1e-3) expressions
        --normalize [bool]   Normalize all expressions
        --fdr [real]     False discovery rate (default: 0.01)  
        --majority_matcher [start|end|complete]     Matching criteria for bindingsite majority vote (only usefull if miRanda, TarPmiR and PITA is enabled)    
      --spongEffects [bool] Use SPONGEs spongEffects extension to identify strong circRNA sponges
        --train [real]  Percentage of samples to use for training (rest will be used for testing) (default = 0.8)
   """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(', ')}"
}

include { CIRCRNA } from "./subworkflows/circRNA.nf"
 
species = params.species ?: params.genome ? params.genomes[ params.genome ].species ?: false : false
fasta = params.fasta ?: params.genome ? params.genomes[ params.genome ].fasta ?: false : false
gtf = params.gtf ?: params.genome ? params.genomes[ params.genome ].gtf ?: false : false
bed12 = params.bed12 ?: params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
miRNA_fasta = params.miRNA_fasta ?: params.genome ? params.genomes[ params.genome ].mature ?: false : false
star = params.star ?: params.genome ? params.genomes[ params.genome ].star ?: false : false

ch_genome_fasta = Channel.fromPath(fasta)
ch_star_index = Channel.fromPath(star)
ch_gtf = Channel.fromPath(gtf)
ch_bed12 = Channel.fromPath(bed12)

workflow CIRCRNA_SPONGING {
    main:
        CIRCRNA(
          ch_star_index,
          ch_gtf,
          ch_genome_fasta,
          ch_bed12
        )
}

workflow {
    CIRCRNA_SPONGING()
}