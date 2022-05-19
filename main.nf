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

nextflow.enable.dsl=1

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

def get_circRNA_paths(LinkedHashMap row) {
    def array = []
    if (!file(row.totalRNA1).exists()) {
        exit 1, "Error: Fastq file does not exist!\n${row.totalRNA1}"
    }
    if (params.single_end) {
        array = [ row.sample, [ file(row.totalRNA1) ] ]
    } else {
        if (!file(row.totalRNA2).exists()) {
             exit 1, "Error: Fastq file does not exist!\n${row.totalRNA2}"
        }
        array = [ row.sample, [ file(row.totalRNA1), file(row.totalRNA2) ] ]
    }
    if (params.differential_expression){
        if (!row.containsKey("condition")) {
            exit 1, "Error: Condition marker missing!"
        }
    }
    return array
}

def get_miRNA_paths(LinkedHashMap row) {
    def array = []
    if (!file(row.smallRNA).exists()) {
        exit 1, "Error: FastQ file does not exist!\n${row.smallRNA}"
    }
    array = [ row.sample, [ file(row.smallRNA) ] ]
    return array
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

// fill params with iGenomes
params.STAR_index = params.STAR_index ?: params.genome ? params.genomes[ params.genome ].star ?: false : false
params.species = params.species ?: params.genome ? params.genomes[ params.genome ].species ?: false : false
params.fasta = params.fasta ?: params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.gtf = params.gtf ?: params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bed12 = params.bed12 ?: params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
params.miRNA_fasta = params.miRNA_fasta ?: params.genome ? params.genomes[ params.genome ].mature ?: false : false
if(!params.miRNA_raw_counts) {
    params.miRNA_related_fasta = params.miRNA_related_fasta ?: params.genome ? params.genomes[ params.genome ].mature_rel ?: false : false
    params.hairpin_fasta = params.hairpin_fasta ?: params.genome ? params.genomes[ params.genome ].hairpin ?: false : false
}

log.info params.species

// create files
Channel.value(file(params.fasta)).into { ch_fasta; ch_fasta_star }
ch_gtf = Channel.value(file(params.gtf))
ch_bed12 = Channel.value(file(params.bed12))
Channel.value(file(params.miRNA_fasta)).into { mirna_fasta_miRanda; mirna_fasta_PITA; mirna_fasta_TarPmiR; mirna_fasta_miRDeep2; mirna_fasta_SPONGE }
ch_miRNA_related_fasta = Channel.value(file(params.miRNA_related_fasta))
ch_hairpin_fasta = Channel.value(file(params.hairpin_fasta))

// Sequencing presets
if (params.protocol == "illumina"){
    params.miRNA_adapter = "TGGAATTCTCGGGTGCCAAGG"
} else if (params.protocol == "nextflex"){
    params.miRNA_adapter = "TGGAATTCTCGGGTGCCAAGG"
} else if (params.protocol == "qiaseq"){
    params.miRNA_adapter = "AACTGTAGGCACCATCAAT"
} else if (params.protocol == "cats"){
    params.miRNA_adapter = "AAAAAAAA"
}

/*
 * CREATE CHANNELS FOR INPUT READ FILES
 */
Channel.fromPath(params.samplesheet)
   .splitCsv( header:true, sep:'\t')
   .map { get_circRNA_paths(it) }.into { ch_totalRNA_reads1; ch_totalRNA_reads2 }

/*
* GENERATE STAR INDEX IN CASE IT IS NOT ALREADY PROVIDED BY iGenomes
*/
process generate_star_index{
    label 'process_high'
    publishDir "${params.outdir}/", mode: params.publish_dir_mode

    input:
    file(fasta) from ch_fasta_star
    file(gtf) from ch_gtf           
    
    output:
    file("star_index") into generated_star_index
                      
    when: (params.STAR_index == null)

    script:
    """
    echo "star index is running"
    mkdir star_index
      
    STAR \\
    --runMode genomeGenerate \\
    --runThreadN 8 \\
    --sjdbGTFfile $gtf \\
    --genomeDir star_index/ \\
    --genomeFastaFiles $fasta
    """
}

ch_star_index = params.STAR_index ? Channel.value(file(params.STAR_index)) : generated_star_index

/*
* PERFORM READ MAPPING OF totalRNA SAMPLES USING STAR
*/
process STAR {
    label 'process_high'
    publishDir "${params.outdir}/samples/${sampleID}/circRNA_detection/", mode: params.publish_dir_mode
    
    input:
    set val(sampleID), file(reads) from ch_totalRNA_reads1
    file star_index from ch_star_index

    output:
    tuple val(sampleID), file("Chimeric.out.junction") into chimeric_junction_files

    script:
    """
    STAR --chimSegmentMin 10 --runThreadN 10 --genomeDir $star_index --readFilesCommand zcat --readFilesIn $reads
    gzip Aligned.out.sam
    """
}

/*
* PARSE STAR OUTPUT INTO CIRCExplorer2 FORMAT
*/
process circExplorer2_Parse {
    label 'process_medium'

    publishDir "${params.outdir}/samples/${sampleID}/circRNA_detection/circExplorer2", mode: params.publish_dir_mode
    
    input:
    tuple val(sampleID), file(chimeric_junction) from chimeric_junction_files

    output:
    tuple val(sampleID), file("back_spliced_junction.bed") into backspliced_junction_bed_files

    script:
    """
    CIRCexplorer2 parse -b "back_spliced_junction.bed" -t STAR $chimeric_junction        
    """
}

/*
* PERFORM circRNA QUANTIFICATION USING CIRCExplorer2
*/
process circExplorer2_Annotate {
    label 'process_medium'

    publishDir "${params.outdir}/samples/${sampleID}/circRNA_detection/circExplorer2", mode: params.publish_dir_mode
    
    input:
    tuple val(sampleID), file(backspliced_junction_bed) from backspliced_junction_bed_files
    file(fasta) from ch_fasta
    file(bed12) from ch_bed12

    output:
    file("${sampleID}_circularRNA_known.txt") into ch_circRNA_known_files

    script:
    """
    CIRCexplorer2 annotate -r bed12 -g $fasta -b $backspliced_junction_bed -o "${sampleID}_circularRNA_known.txt"
    """
}

/*
* MERGE RAW circRNA RESULTS INTO ONE TABLE SUMMARIZING ALL SAMPLES
*/
process summarize_detected_circRNAs{
    label 'process_medium'

    publishDir "${params.outdir}/results/circRNA/", mode: params.publish_dir_mode
    
    input:
    file(circRNA_file) from ch_circRNA_known_files.collect()

    output:
    file("circRNA_counts_raw.tsv") into (ch_circRNA_counts_raw1, ch_circRNA_counts_raw2)

    script:
    """
    Rscript "${projectDir}"/bin/circRNA_summarize_results.R $params.samplesheet $params.outdir
    """
}

/*
* FOR THE PREVIOUSLY DETECTED circRNAs EXTRACT FASTA SEQUENCES ACOORDING TO circRNA TYPE
*/
process extract_circRNA_sequences {
    label 'process_medium'
    publishDir "${params.outdir}/results/binding_sites/input/", mode: params.publish_dir_mode
    
    input:
    file(circRNAs_raw) from ch_circRNA_counts_raw1

    output:
    file("circRNAs.fa") into circRNAs_raw_fasta

    script:
    """
	Rscript "${projectDir}"/bin/extract_fasta.R $params.genome $circRNAs_raw
    """
}

/*
* CREATE PSIRC INDEX IF NOT ALREADY PRESENT/ GIVEN
*/
psirc = params.psirc_exc ? params.psirc_exc : "psirc-quant"
psirc_index_path = params.psirc_index ? params.psirc_index : params.outdir + "/results/psirc/psirc.index"
if(!file(psirc_index_path).exists()) {
    process psirc_index {
        label 'process_medium'
        publishDir "${params.outdir}/results/psirc/", mode: params.publish_dir_mode

        input:
        file(circ_fasta) from circRNAs_raw_fasta
        val(psirc_quant) from psirc

        output:
        file("psirc.index") into psirc_index

        script:
        """
        Rscript "${projectDir}"/bin/build_psirc_index.R \
        --index "psirc.index" \
        --transcriptome $params.transcriptome \
        --circ_fasta $circ_fasta \
        --psirc_quant $psirc_quant
        """
    }
} else {
    Channel.value( file(psirc_index_path) ).into{ psirc_index }
}

/*
* QUANTIFY EXPRESSIONS USING PSIRC
*/
psirc_out = params.outdir + "/results/psirc/"
if (!file(psirc_out + "quant_linear_expression.tsv").exists()) {
    process psirc_quant {
        label 'process_medium'
        publishDir "${params.outdir}/results/psirc/tmp/", mode: params.publish_dir_mode

        input:
        set val(sampleID), file(reads) from ch_totalRNA_reads2
        val(psirc_quant) from psirc
        val(psirc_index) from psirc_index

        output:
        file("${sampleID}/abundance.tsv") into psirc_outputs

        script:
        if (params.single_end)
            """
            $psirc_quant quant -i $psirc_index -o $sampleID --single -l 76 -s 20 $reads
            """
        else
            """
            $psirc_quant quant -i $psirc_index -o $sampleID $reads
            """
    }

    /*
    * COMBINE PSIRC OUTPUTS -> quantified linear and circular expression for each sample
    */
    process process_psirc {
        label 'process_medium'
        publishDir "${params.outdir}/results/psirc/", mode: params.publish_dir_mode

        input:
        file(circ_counts) from ch_circRNA_counts_raw2
        file("abundance.tsv") from psirc_outputs.collect()

        output:
        file("quant_circ_expression.tsv") into ch_circRNA_counts_raw_quant
        file("quant_linear_expression.tsv") into (gene_expression1, gene_expression2)
        file("TPM_map.tsv") into (TPM_map1, TPM_map2)
        file("quant_effects.png") into quant_effects

        script:
        """
        Rscript "${projectDir}"/bin/quantify_circ_expression.R \
        --circ_counts $circ_counts \
        --dir "${params.outdir}/results/psirc/tmp/" \
        --samplesheet $params.samplesheet \
        --pseudocount $params.pseudocount
        """
    }
} else {
    Channel.fromPath(psirc_out + "quant_circ_expression.tsv").into{ ch_circRNA_counts_raw_quant }
    Channel.fromPath(psirc_out + "quant_linear_expression.tsv").into{ gene_expression1; gene_expression2 }
    Channel.fromPath(psirc_out + "TPM_map.tsv").into{ TPM_map1; TPM_map2 }
}
// choose either quantified or regular circRNA counts for downstream analysis
if (params.quantification){
    ch_circRNA_counts_raw_quant.into{ ch_circRNA_counts1; ch_circRNA_counts2 }
} else {
    ch_circRNA_counts_raw2.into{ ch_circRNA_counts1; ch_circRNA_counts2 }
}

/*
* NORMALIZE RAW circRNA COUNT USING LIBRARY SIZE ESTIMATION
*/
process normalize_circRNAs{
    label 'process_medium'

    publishDir "${params.outdir}/results/circRNA/", mode: params.publish_dir_mode
    
    input:
    file(circRNA_counts_raw) from ch_circRNA_counts1

    output:
    file("circRNA_counts_normalized.tsv") into (ch_circRNA_counts_norm1, ch_circRNA_counts_norm2)

    script:
    """
    Rscript "${projectDir}"/bin/circRNA_results_LibrarySizeEstimation.R $circRNA_counts_raw $params.samplesheet $params.outdir
    """
}

/*
* FILTER circRNAs TO REDUCE LOW EXPRESSED ONES
*/
process filter_circRNAs{
    label 'process_medium'

    publishDir "${params.outdir}/results/circRNA/", mode: params.publish_dir_mode
    
    input:
    file(circRNA_counts_norm) from ch_circRNA_counts_norm1

    output:
    file("circRNA_counts_filtered.tsv") into ch_circRNA_counts_filtered

    script:
    """
    Rscript "${projectDir}"/bin/circRNA_filtering.R $circRNA_counts_norm $params.samplesheet $params.outdir $params.sample_percentage $params.read_threshold
    """
}

/*
* DATABASE ANNOTATION USING LIFTOVER FOR GENOMIC COORDINATE CONVERSION AND CIRCBASE
*/
if (params.database_annotation){
    circ_annotation = params.outdir + "/results/circRNA/circRNAs_annotated.tsv"
    circ_counts_annotated_path = params.outdir + "/results/circRNA/circRNA_counts_annotated.tsv"
    if (!file(circ_counts_annotated_path).exists()) {
        process database_annotation{
        label 'process_medium'

        publishDir "${params.outdir}/results/circRNA/", mode: params.publish_dir_mode

        input:
        file(circRNAs_filtered) from ch_circRNA_counts_filtered

        output:
        file("circRNAs_annotated.tsv") into circRNAs_annotated
        file("circRNA_counts_annotated.tsv") into (ch_circRNA_counts_filtered1, ch_circRNA_counts_filtered2, ch_circRNA_counts_filtered3, ch_circRNA_counts_filtered4, ch_circRNA_counts_filtered5)
        script:
        if( params.offline_circ_db == null )
            """
            python3 "${projectDir}"/bin/circRNA_db_annotation.py -o $params.species -gv $params.genome -d $circRNAs_filtered -ao $params.annotated_only
            """
        else
            """
            python3 "${projectDir}"/bin/circRNA_db_annotation.py -o $params.species -gv $params.genome -d $circRNAs_filtered -off $params.offline_circ_db -ao $params.annotated_only
            """
        }
    } else {
        Channel.fromPath(circ_annotation).into{ circRNAs_annotated }
        Channel.fromPath(circ_counts_annotated_path).into{ ch_circRNA_counts_filtered1; ch_circRNA_counts_filtered2; ch_circRNA_counts_filtered3; ch_circRNA_counts_filtered4; ch_circRNA_counts_filtered5 }
    }
} else {
    ch_circRNA_counts_filtered.into{ ch_circRNA_counts_filtered1; ch_circRNA_counts_filtered2; ch_circRNA_counts_filtered3; ch_circRNA_counts_filtered4; ch_circRNA_counts_filtered5 }
}

/*
* EXTRACT circRNA FASTAS
*/
process circ_fastas{
    label 'process_medium'

    publishDir "${params.outdir}/results/binding_sites/input/", mode: params.publish_dir_mode

    input:
    file(circRNAs_filtered) from ch_circRNA_counts_filtered1

    output:
    file("circRNAs.fa") into (circRNAs_fasta1, circRNAs_fasta2, circRNAs_fasta3)

    script:
    """
	Rscript "${projectDir}"/bin/extract_fasta.R $params.genome $circRNAs_filtered
    """
}

/*
* DIFFERENTIAL EXPRESSION ANALYSIS USING SAM FILES FROM STAR
*/
if (params.differential_expression){
    process differential_expression {
        label 'process_medium'
        publishDir "${params.outdir}/results/differential_expression", mode: params.publish_dir_mode
        errorStrategy 'ignore'

        input:
        file(circRNA_counts) from ch_circRNA_counts_filtered2
        file(circRNA_raw) from ch_circRNA_counts2
        file(gene_expression) from gene_expression1
        val(tpm_map) from TPM_map1

        output:
        file("total_rna/total_rna.tsv") into deseq_total_rna
        file("total_rna/total_rna.signif.tsv") into DE_mRNA_signif
        file("circ_rna_DE/circ_rna_DE.tsv") into deseq_circ_rna
        file("circ_rna_DE/circ_rna_DE.signif.tsv") into DE_circ_signif
        file("total_rna/*.png") into total_plots
        file("circ_rna_DE/*.png") into circ_plots
        file("*.png") into supplementary_plots
        file("DESeq2.RData") into deseq2_rdata

        script:
        """
        Rscript "${projectDir}"/bin/differentialExpression.R \\
            --gene_expr $gene_expression \\
            --samplesheet $params.samplesheet \\
            --circ_filtered $circRNA_counts \\
            --circ_raw $circRNA_raw \\
            --tpm_map $tpm_map
        """
    }
}

/*
* DETERMINE miRNA BINDING SITES ON THE PREVIOUSLY DETECTED circRNAs USING miranda
run only if file is not already present
*/
miranda_path = params.outdir + "/results/binding_sites/output/miRanda"
miranda_output = miranda_path + "/bind_sites_raw.out"
if (!file(miranda_output).exists()) {
    miranda_tmp = miranda_path + "/tmp"
    process miRanda {
        label 'process_medium'
        publishDir miranda_tmp, mode: params.publish_dir_mode
        
        input:
        file(circRNA_fasta) from circRNAs_fasta1.splitFasta( by: params.splitter, file: true )
        file(miRNA_fasta) from mirna_fasta_miRanda

        output:
        file("bind_sites_raw.out") into bind_sites_split

        script:
        """
        miranda $miRNA_fasta $circRNA_fasta -out "bind_sites_raw.out" -quiet
        """
    }
    // combine files to one
    bind_sites_split.collectFile(name: miranda_output, newLine: true).into{ bind_sites_out }
    // delete tmp files
    file(miranda_tmp).deleteDir()
} else {
    Channel.fromPath(miranda_output).into{ bind_sites_out }
}
/*
* PROCESS miranda OUTPUT INTO A TABLE FORMAT
*/
process binding_sites_processing {
    label 'process_medium'
    publishDir miranda_path, mode: params.publish_dir_mode
    
    input:
    file(bind_sites_raw) from bind_sites_out

    output:
    file("bind_sites_processed.txt") into bind_sites_processed

    script:
    """
    echo -e "miRNA\tTarget\tScore\tEnergy-Kcal/Mol\tQuery-Al(Start-End)\tSubject-Al(Start-End)\tAl-Len\tSubject-Identity\tQuery-Identity" > "bind_sites_processed.txt"
    grep -A 1 "Scores for this hit:" $bind_sites_raw | sort | grep ">" | cut -c 2- >> "bind_sites_processed.txt"
    """
}

/*
* FILTER BINDING SITES, KEEP TOP 25%
*/
process binding_sites_filtering {
    label 'process_medium'
    publishDir miranda_path, mode: params.publish_dir_mode
    
    input:
    file(bind_sites_proc) from bind_sites_processed
    
    output:
    file("bindsites_25%_filtered.tsv") into (ch_bindsites_filtered1, ch_bindsites_filtered2)

    script:
    """
    Rscript "${projectDir}"/bin/binding_sites_analysis.R ${bind_sites_proc}
    """
}

/*
* RUN TARPMIR ANALYSIS ON CIRCRNA FASTAS
*/
if (params.tarpmir) {
    model = params.model ? Channel.value(file(params.model)) : Channel.value(file(projectDir + "/data/tarpmir_models/Human_sklearn_0.22.pkl"))
    tarpmir_path = "${params.outdir}/results/binding_sites/output/TarPmiR"
    tarpmir_tmp = tarpmir_path + "/tmp"
    // RUN TARPMIR ON CHUNKED MRNA FASTAS
    tarpmir_out = tarpmir_path + "/TarPmiR_bp.tsv"
    if (!file(tarpmir_out).exists()) {
        process TarPmiR {
            label 'process_medium'
            publishDir tarpmir_tmp, mode: params.publish_dir_mode

            input:
            file(model) from model
            file(mRNA_fasta) from circRNAs_fasta2.splitFasta( by: params.splitter, file: true )
            file(miRNA_fasta) from mirna_fasta_TarPmiR

            output:
            file("bindings.bp") into bp_files

            script:
            """
            python3 "${projectDir}/bin/TarPmiR_threading.py" \\
            -a $miRNA_fasta \\
            -b $mRNA_fasta \\
            -m $model \\
            -p $params.p \\
            -t $params.threads \\
            -o "bindings.bp"
            """
        }

        // combine files to one
        bp_files.collectFile(name: tarpmir_out, newLine: true).into{ tarpmir_bp_file }
        // delete tmp files
        file(tarpmir_tmp).deleteDir()
    } else {
        Channel.fromPath(tarpmir_out).into{ tarpmir_bp_file }
    }
} else {
    Channel.of( 'null' ).into{ tarpmir_bp_file }
}

/*
 * RUN PITA ANALYSIS FOR circRNAs
 */
if (params.pita) {
    pita_tmp = "${params.outdir}/results/binding_sites/output/PITA/tmp"
    pita_out = "${params.outdir}/results/binding_sites/output/PITA/circRNA_pita_results.tsv"
    pita_options = "-l " + params.pita_l + " -gu " + params.pita_gu + " -m " + params.pita_m
    if(!file(pita_out).exists()){
        process PITA {
            label 'process_long'
            publishDir pita_tmp, mode: params.publish_dir_mode
            errorStrategy 'retry'

            input:
            file(circ_fasta) from circRNAs_fasta3.splitFasta( by: params.splitter, file: true )
            val(options) from pita_options
            file(miRNA_fasta) from mirna_fasta_PITA

            output:
            file("circRNA_pita_results.tab") into pita_splits

            script:
            """
            pita_prediction.pl -utr $circ_fasta -mir $miRNA_fasta -prefix circRNA $options
            """
        }

        // collect all PITA splits
        pita_splits.collectFile(name: pita_out, newLine: true).into{ pita_results }
        // delete tmp directory
        file(pita_tmp).deleteDir()
    } else {
        Channel.fromPath(pita_out).into{ pita_results }
    }
} else {
    Channel.of( 'null' ).into{ pita_results }
}

/*
 * RUN miRNA part only if circRNA_only==false
 */

if (!params.circRNA_only) {
    /*
    * GET miRNA RAW COUNTS
    */
    if( params.miRNA_raw_counts != null ) {

        /*
        * IF RAW miRNA COUNTS ARE ALREADY SPECIFIED IN A FILE
        */
        ch_miRNA_counts_raw = Channel.fromPath(params.miRNA_raw_counts) 

    } else {
    
        /*
        * PERFORM miRNA DETECTION USING miRDeep2 FROM SPECIFIED READ FILES
        * CREATE INPUT CHANNEL
        */

        ch_smallRNA_reads=Channel.fromPath(params.samplesheet)
        .splitCsv( header:true, sep:'\t')
        .map { get_miRNA_paths(it) }

        /*
        * GENERATE BOWTIE INDEX IN CASE IT IS NOT ALREADY PROVIDED
        */
        process generate_bowtie_index{
            label 'process_high'
            publishDir "${params.outdir}/bowtie_index/", mode: params.publish_dir_mode

            input:
            file(fasta) from ch_fasta
            
            output:
            file("${fasta.baseName}*") into ch_generated_bowtie_index_files
            val("${fasta.baseName}") into ch_generated_bowtie_index
                            
            when: (params.bowtie_index == null)

            script:
            """
            echo "bowtie index is in ${fasta.baseName}"
            bowtie-build $fasta ${fasta.baseName}
            """
        }

        ch_bowtie_index = params.bowtie_index ? Channel.value(params.bowtie_index) : ch_generated_bowtie_index

        /*
        * PERFORM miRNA READ MAPPING USING miRDeep2
        */
        miRNA_adapter = params.miRNA_adapter ? "-k " + params.miRNA_adapter : ""
        process miRDeep2_mapping {
            label 'process_high'
            publishDir "${params.outdir}/samples/${sampleID}/miRNA_detection/", mode: params.publish_dir_mode

            input:
            tuple val(sampleID), file(read_file) from ch_smallRNA_reads
            val(index) from ch_bowtie_index
            val(adapter) from miRNA_adapter

            output: 
            tuple val(sampleID), file("reads_collapsed.fa"), file("reads_vs_ref.arf") into ch_miRNA_mapping_output

            script:
            """
            gunzip < $read_file > "${sampleID}.fastq"
            mapper.pl "${sampleID}.fastq" -e -h -i $adapter -l 18 -m -p $index -s "reads_collapsed.fa" -t "reads_vs_ref.arf" -v
            """
        }

        /*
        * PERFORM miRNA QUANTIFICATION USING miRDeep2
        */
        process miRDeep2_quantification {
            label 'process_high'
            publishDir "${params.outdir}/samples/${sampleID}/miRNA_detection/", mode: params.publish_dir_mode
    
            input:
            tuple val(sampleID), file(reads_collapsed_fa), file(reads_vs_ref_arf) from ch_miRNA_mapping_output
            file(fasta) from ch_fasta
            file(miRNA_fasta) from mirna_fasta_miRDeep2
            file(miRNA_related_fasta) from ch_miRNA_related_fasta
            file(hairpin_fasta) from ch_hairpin_fasta

            output:
            file("miRNAs_expressed*") into ch_miRNA_expression_files

            script:
            """
            miRDeep2.pl $reads_collapsed_fa $fasta $reads_vs_ref_arf $miRNA_fasta $miRNA_related_fasta $hairpin_fasta -t $params.species -d -v 
            """
        }

        /*
        * MERGE RAW miRNA RESULTS INTO ONE TABLE SUMMARIZING ALL SAMPLES
        */
        process summarize_detected_miRNAs{
            label 'process_medium'

            publishDir "${params.outdir}/results/miRNA/", mode: params.publish_dir_mode
        
            input:
            val(miRNAs_expressed) from ch_miRNA_expression_files.unique().collect()

            output:
            file("miRNA_counts_raw.tsv") into ch_miRNA_counts_raw

            script:
            """
            Rscript "${projectDir}"/bin/miRNA_summarize_results.R $params.samplesheet $params.outdir
            """
        }
    }

    /*
    * NORMALIZE RAW miRNA COUNT USING LIBRARY SIZE ESTIMATION
    */
    if (params.miRNA_normalization){
        process normalize_miRNAs{
            label 'process_low'

            publishDir "${params.outdir}/results/miRNA/", mode: params.publish_dir_mode
            
            input:
            file(miRNA_counts_raw) from ch_miRNA_counts_raw

            output:
            file("miRNA_counts_normalized.tsv") into (ch_miRNA_counts_norm1, ch_miRNA_counts_norm2)

            script:
            """
            Rscript "${projectDir}"/bin/miRNA_results_LibrarySizeEstimation.R $miRNA_counts_raw $params.outdir
            """
        }
    } else {
        ch_miRNA_counts_raw.into{ ch_miRNA_counts_norm1; ch_miRNA_counts_norm2 }
    }
    

    /*
    * FILTER miRNAs TO REDUCE LOW EXPRESSED ONES
    */
    if (params.miRNA_filtering) {
        process filter_miRNAs{
            label 'process_medium'

            publishDir "${params.outdir}/results/miRNA/", mode: params.publish_dir_mode
            
            input:
            file(miRNA_counts_norm) from ch_miRNA_counts_norm1

            output:
            file("miRNA_counts_filtered.tsv") into (ch_miRNA_counts_filtered1, ch_miRNA_counts_filtered2, ch_miRNA_counts_filtered3)

            script:
            """
            Rscript "${projectDir}"/bin/miRNA_filtering.R $miRNA_counts_norm $params.outdir $params.sample_percentage $params.read_threshold
            """
        }
    } else {
        ch_miRNA_counts_norm1.into{ch_miRNA_counts_filtered1; ch_miRNA_counts_filtered2; ch_miRNA_counts_filtered3}
    }
    
    /*
    * FOR ALL POSSIBLE circRNA-miRNA PAIRS COMPUTE PEARSON CORRELATION
    * USES ONLY miRanda bindsites
    */
    if (params.correlations){
        process compute_correlations{
            label 'process_high'

            publishDir "${params.outdir}/results/sponging/", mode: params.publish_dir_mode
            
            input:
            file(miRNA_counts_filtered) from ch_miRNA_counts_filtered1
            file(circRNA_counts_filtered) from ch_circRNA_counts_filtered3
            file(filtered_bindsites) from ch_bindsites_filtered1

            output:
            file("filtered_circRNA_miRNA_correlation.tsv") into ch_correlations

            script:
            """
            Rscript "${projectDir}"/bin/compute_correlations.R $params.samplesheet $miRNA_counts_filtered $circRNA_counts_filtered $filtered_bindsites
            """
        }

        /*
        * ANALYZE THE CORRELATION OF ALL PAIRS AND DETERMINE OVERALL DISTRIBUTION
        * USING BINDING SITES INFORMATION. COMPUTE STATISTICS AND PLOTS
        */
        process correlation_analysis{
            label 'process_high'

            publishDir "${params.outdir}/results/sponging/", mode: params.publish_dir_mode
            
            input:
            file(correlations) from ch_correlations
            file(miRNA_counts_filtered) from ch_miRNA_counts_filtered2
            file(circRNA_counts_filtered) from ch_circRNA_counts_filtered4
            file(miRNA_counts_norm) from ch_miRNA_counts_norm2
            file(circRNA_counts_norm) from ch_circRNA_counts_norm2


            output:
            file("sponging_statistics.txt") into ch_sponging_statistics
            file("plots/*.png") into ch_plots

            script:
            """
            mkdir -p "${params.outdir}/results/sponging/plots/"
            Rscript "${projectDir}"/bin/correlation_analysis.R $params.samplesheet $miRNA_counts_filtered $circRNA_counts_filtered $correlations $params.outdir $params.sample_group $miRNA_counts_norm $circRNA_counts_norm
            """
        }
    }

    /*
    * SPONGE ANALYSIS (https://github.com/biomedbigdata/SPONGE)
    */
    if (params.sponge) {
        // USE GIVEN TARGET SYMBOLS OR DEFAULT LOCATED IN DATA
        target_scan_symbols = params.target_scan_symbols ? Channel.value(file(params.target_scan_symbols)) : Channel.value(file(projectDir + "/data/miRNA_target_symbols/hsa_mirWalk_lncbase_21_ENSG.tsv.gz"))
        process SPONGE {
            label 'process_long'
            errorStrategy 'ignore'

            publishDir "${params.outdir}/results/sponging/SPONGE", mode: params.publish_dir_mode

            input:
            file(gene_expression) from gene_expression2
            file(circRNA_counts_filtered) from ch_circRNA_counts_filtered5
            file(mirna_expression) from ch_miRNA_counts_filtered3
            file(miRNA_fasta) from mirna_fasta_SPONGE
            file(miranda_bind_sites) from ch_bindsites_filtered2
            file(target_scan_symbols) from target_scan_symbols
            val(tarpmir) from tarpmir_bp_file
            val(pita) from pita_results
            val(normalize) from params.normalize ? "--normalize" : ""
            val(tpm_map) from TPM_map2
            val(tpm) from params.tpm ? "--tpm" : ""

            output:
            file("sponge.RData") into (sponge_rimage1, sponge_rimage2)
            file("plots/*.png") into plots
            file("circRNA/*") into circResults
            file("total/*") into totalResults

            script:
            """
            Rscript "${projectDir}"/bin/SPONGE.R \\
            --gene_expr $gene_expression \\
            --circ_rna $circRNA_counts_filtered \\
            --mirna_expr $mirna_expression \\
            --mir_fasta $miRNA_fasta \\
            --fdr $params.fdr \\
            --target_scan_symbols $target_scan_symbols \\
            --miranda_data $miranda_bind_sites \\
            --tarpmir_data $tarpmir \\
            --pita_data $pita \\
            --majority_matcher $params.majority_matcher \\
            --pseudocount $params.pseudocount \\
            --tpm_map $tpm_map \\
            $tpm \\
            $normalize
            """
        }
        /*
        * PERFORM SPONGE circRNA DE ANALSIS
        */
        if (params.differential_expression) {
            process SPONGE_DE_ANALYSIS {
                label 'process_low'
                errorStrategy 'ignore'

                publishDir "${params.outdir}/results/sponging/SPONGE/DE_analysis", mode: params.publish_dir_mode

                input:
                file(sponge_data) from sponge_rimage1
                file(circ_signif_DE) from DE_circ_signif
                file(mRNA_signif_DE) from DE_mRNA_signif

                output:
                file("DE_SPONGE.html") into DE_SPONGE_graph

                script:
                """
                Rscript "${projectDir}"/bin/SPONGE_analysis.R \\
                $sponge_data \\
                $circ_signif_DE \\
                $mRNA_signif_DE
                """
            }
        }
        /*
        * RUN spongEffects ON SPONGE RESULTS
        */
        if (params.spongEffects) {
            process spongEffects {
                label 'process_high'
                errorStrategy 'ignore'

                publishDir "${params.outdir}/results/sponging/SPONGE/spongEffects", mode: params.publish_dir_mode

                input:
                file(sponge_data) from sponge_rimage2

                output:
                file("*") into spongEffectsResults

                script:
                """
                Rscript "${projectDir}"/bin/spongEffects.R \\
                --spongeData $sponge_data \\
                --meta $params.samplesheet \\
                --train $params.se_train
                """
            }
        }
    }
}
