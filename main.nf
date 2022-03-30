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

    nextflow run nf-core/circrnasponging 

    Mandatory arguments:
      --samplesheet [file]		Path to samplesheet (must be surrounded with quotes)
      --out_dir [file]			The output directory where the results will be saved
      --species [str]			Species name in 3 letter code (hsa for human, mmu for mouse)
      --genome_version [str]    Genome version that will be used for mapping e.g. for human hg19 or hg38
      --miRNA_adapter [str] 		miRNA adapter used for trimming
      -profile [str]           	 	Configuration profile to use. Can use multiple (comma separated)
                                      	Available: conda, docker, singularity, test, awsbatch, <institute> and more

      --fasta [file] 			Path to genome fasta (must be surrounded with quotes)
      --gtf [file]			Path to gtf file (must be surrounded with quotes)
      --gene_pred [file]		Path to gene annotation (must be surrounded with quotes)
      --mature_fasta [file]		Path to mature miRNA fasta (must be surrounded with quotes)
      --mature_other_fasta [file]	Path to mature miRNA fasta of related species (must be surrounded with quotes)
      --hairpin_fasta [file]		Path to miRNA hairpin fasta (must be surrounded with quotes)
      if SPONGE is enabled:
        Supply at least one of the following target scan symbols:
            --target_scan_symbols [file]    Path to target scan symbols in tsv and SPONGE format (rows=GENE,cols=miRNA,data=counts)
            --lncBaseData    [file]     Path to lncBase targets in tsv
            --miRTarBaseData [file]     Path to miRTarBase targets in tsv
            --TargetScanData [file]     Path to TargetScan targets in tsv
            --miRDB_data     [file]     Path to miRDB targets in tsv
        --fdr   [real]     False discovery rate (default: 0.1)

    Options:
      --miRNA_raw_counts [file]		Path to tabulated raw miRNA counts (must be surrounded with quotes)
      --single_end [bool]            	Specifies that the total RNA input is single-end reads
      --sample_group	[file]		File specifying partitioning of samples into groups (must be surrounded with quotes)
      --read_threshold [real]		Positive. Read counts under this threshold are considered to be low expressed
      --sample_percentage [real]	Between 0 and 1. Minimum percentage of samples that should have no low expression
      --circRNA_only [bool]  		Run only circRNA analysis, don't run miRNA analysis
      --database_annotation [bool]  Annotate circRNA hits with circBase data
        --offline_circ_db [file]      File containing downloaded circBase entries for offline access to the database
      --differential_expression [bool]  Enable differential expression analysis using DESeq2 on all given RNA-seq data and circRNA only
      --tarpmir [bool]   Wheather to use tarpmir for bindsite prediction (may take significantly longer)
        --model [pkl]       Path to a specific tarpmir model as pickle (default model is Human_sklearn_0.22.pkl, located in data directory)
        --p     [real]      Double between 0 and 1 specifing cutoff for tarpmir
        --threads [int]     Number of threads to use for tarpmir
        --splitter [int]    Number of fasta circRNA entries to use for each tarpmir run
      --sponge [bool]   Wheather to perform SPONGE ceRNA network analysis          
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

def check_input(){
    if (!params.genome_version && params.database_annotation) {
        exit 1, "Error: genome version not specified, which is mandatory for database annotation"
    }
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * CREATE CHANNELS FOR INPUT READ FILES
 */
Channel.fromPath(params.samplesheet)
   .splitCsv( header:true, sep:'\t')
   .map { get_circRNA_paths(it) }.into { ch_totalRNA_reads1; ch_totalRNA_reads2 }

ch_fasta = Channel.value(file(params.fasta))
ch_gtf = Channel.value(file(params.gtf))

/*
* CHECK INPUT OPTIONS
*/
check_input()

/*
* GENERATE STAR INDEX IN CASE IT IS NOT ALREADY PROVIDED
*/
process generate_star_index{
    label 'process_high'
    publishDir "${params.out_dir}/", mode: params.publish_dir_mode

    input:
    file(fasta) from ch_fasta
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

/*
* GENERATE SALMON INDEX FOR GIVEN ORGANISM
*/
process generate_salmon_index {
    label 'process_high'
    publishDir "${params.out_dir}/", mode: params.publish_dir_mode

    output:
    file("salmon_index") into generated_salmon_index

    when: (params.SALMON_index == null && params.transcriptome != null)

    script:
    """
    salmon index -t $params.transcriptome -i salmon_index
    """
}

ch_star_index = params.STAR_index ? Channel.value(file(params.STAR_index)) : generated_star_index

ch_salmon_index = params.SALMON_index ? Channel.value(file(params.SALMON_index)) : generated_salmon_index

/*
* PERFORM READ MAPPING OF totalRNA SAMPLES USING STAR
*/
process STAR {
    label 'process_high'
    publishDir "${params.out_dir}/samples/${sampleID}/circRNA_detection/", mode: params.publish_dir_mode
    
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

    publishDir "${params.out_dir}/samples/${sampleID}/circRNA_detection/circExplorer2", mode: params.publish_dir_mode
    
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

    publishDir "${params.out_dir}/samples/${sampleID}/circRNA_detection/circExplorer2", mode: params.publish_dir_mode
    
    input:
    tuple val(sampleID), file(backspliced_junction_bed) from backspliced_junction_bed_files
    file(fasta) from ch_fasta

    output:
    file("${sampleID}_circularRNA_known.txt") into ch_circRNA_known_files

    script:
    """
    CIRCexplorer2 annotate -r $params.gene_pred -g $fasta -b $backspliced_junction_bed -o "${sampleID}_circularRNA_known.txt"
    """
}

/*
* MERGE RAW circRNA RESULTS INTO ONE TABLE SUMMARIZING ALL SAMPLES
*/
process summarize_detected_circRNAs{
    label 'process_medium'

    publishDir "${params.out_dir}/results/circRNA/", mode: params.publish_dir_mode
    
    input:
    file(circRNA_file) from ch_circRNA_known_files.collect()

    output:
    file("circRNA_counts_raw.tsv") into (ch_circRNA_counts_raw1, ch_circRNA_counts_raw2)

    script:
    """
    Rscript "${projectDir}"/bin/circRNA_summarize_results.R $params.samplesheet $params.out_dir
    """
}

/*
* FOR THE PREVIOUSLY DETECTED circRNAs EXTRACT FASTA SEQUENCES ACOORDING TO circRNA TYPE
*/
process extract_circRNA_sequences {
    label 'process_medium'
    publishDir "${params.out_dir}/results/binding_sites/input/", mode: params.publish_dir_mode
    
    input:
    file(circRNAs_raw) from ch_circRNA_counts_raw1

    output:
    file("circRNAs.fa") into circRNAs_raw_fasta

    script:
    """
	Rscript "${projectDir}"/bin/extract_fasta.R $params.genome_version $circRNAs_raw
    """
}

/*
* CREATE PSIRC INDEX IF NOT ALREADY PRESENT/ GIVEN
*/
psirc = params.psirc_exc ? params.psirc_exc : projectDir + "/ext/psirc-quant"
psirc_index_path = params.psirc_index ? params.psirc_index : params.out_dir + "/results/psirc/psirc.index"
if(!file(psirc_index_path).exists()) {
    process psirc_index {
        label 'process_medium'
        publishDir "${params.out_dir}/results/psirc/", mode: params.publish_dir_mode

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
process psirc_quant {
    label 'process_medium'
    publishDir "${params.out_dir}/results/psirc/tmp/", mode: params.publish_dir_mode

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
    publishDir "${params.out_dir}/results/psirc/", mode: params.publish_dir_mode

    input:
    file(circ_counts) from ch_circRNA_counts_raw2
    file("abundance.tsv") from psirc_outputs.collect()

    output:
    file("quant_circ_expression.tsv") into ch_circRNA_counts_raw_quant
    file("quant_linear_expression.tsv") into gene_expression
    file("TPM_map.tsv") into TPM_map
    file("quant_effects.png") into quant_effects

    script:
    """
    Rscript "${projectDir}"/bin/quantify_circ_expression.R \
    --circ_counts $circ_counts \
    --dir "${params.out_dir}/results/psirc/tmp/" \
    --samplesheet $params.samplesheet
    """
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

    publishDir "${params.out_dir}/results/circRNA/", mode: params.publish_dir_mode
    
    input:
    file(circRNA_counts_raw) from ch_circRNA_counts1

    output:
    file("circRNA_counts_normalized.tsv") into (ch_circRNA_counts_norm1, ch_circRNA_counts_norm2)

    script:
    """
    Rscript "${projectDir}"/bin/circRNA_results_LibrarySizeEstimation.R $circRNA_counts_raw $params.samplesheet $params.out_dir
    """
}

/*
* FILTER circRNAs TO REDUCE LOW EXPRESSED ONES
*/
process filter_circRNAs{
    label 'process_medium'

    publishDir "${params.out_dir}/results/circRNA/", mode: params.publish_dir_mode
    
    input:
    file(circRNA_counts_norm) from ch_circRNA_counts_norm1

    output:
    file("circRNA_counts_filtered.tsv") into ch_circRNA_counts_filtered

    script:
    """
    Rscript "${projectDir}"/bin/circRNA_filtering.R $circRNA_counts_norm $params.samplesheet $params.out_dir $params.sample_percentage $params.read_threshold
    """
}

/*
* DATABASE ANNOTATION USING LIFTOVER FOR GENOMIC COORDINATE CONVERSION AND CIRCBASE
*/
if (params.database_annotation){
    process database_annotation{
    label 'process_medium'

    publishDir "${params.out_dir}/results/circRNA/", mode: params.publish_dir_mode

    input:
    file(circRNAs_filtered) from ch_circRNA_counts_filtered

    output:
    file("circRNAs_annotated.tsv") into circRNAs_annotated
    file("circRNA_counts_annotated.tsv") into (ch_circRNA_counts_filtered1, ch_circRNA_counts_filtered2, ch_circRNA_counts_filtered3, ch_circRNA_counts_filtered4, ch_circRNA_counts_filtered5)
    script:
    if( params.offline_circ_db == null )
        """
        python3 "${projectDir}"/bin/circRNA_db_annotation.py -o $params.species -gv $params.genome_version -d $circRNAs_filtered -ao $params.annotated_only
        """
    else
        """
        python3 "${projectDir}"/bin/circRNA_db_annotation.py -o $params.species -gv $params.genome_version -d $circRNAs_filtered -off $params.offline_circ_db -ao $params.annotated_only
        """
    }
} else {
    ch_circRNA_counts_filtered.into{ ch_circRNA_counts_filtered1; ch_circRNA_counts_filtered2; ch_circRNA_counts_filtered3; ch_circRNA_counts_filtered4; ch_circRNA_counts_filtered5 }
}

/*
* EXTRACT circRNA FASTAS
*/
process circ_fastas{
    label 'process_medium'

    publishDir "${params.out_dir}/results/binding_sites/input/", mode: params.publish_dir_mode

    input:
    file(circRNAs_filtered) from ch_circRNA_counts_filtered1

    output:
    file("circRNAs.fa") into (circRNAs_fasta1, circRNAs_fasta2, circRNAs_fasta3)

    script:
    """
	Rscript "${projectDir}"/bin/extract_fasta.R $params.genome_version $circRNAs_filtered
    """
}

/*
* DIFFERENTIAL EXPRESSION ANALYSIS USING SAM FILES FROM STAR
*/
if (params.differential_expression){
    process differential_expression {
        label 'process_medium'
        publishDir "${params.out_dir}/results/gene_expression/differential_expression", mode: params.publish_dir_mode
        errorStrategy 'ignore'

        input:
        file(circRNA_counts) from ch_circRNA_counts_filtered2
        file(circRNA_raw) from ch_circRNA_counts2
        file(gene_expression) from gene_expression

        output:
        file("total_rna/total_rna.tsv") into deseq_total_rna
        file("total_rna/total_rna.signif.tsv") into DE_mRNA_signif
        file("circ_rna_DE/circ_rna_DE.tsv") into deseq_circ_rna
        file("circ_rna_DE/circ_rna_DE.signif.tsv") into DE_circ_signif
        file("total_rna/*.png") into total_plots
        file("circ_rna_DE/*.png") into circ_plots
        file("DESeq2.RData") into deseq2_rdata

        script:
        """
        Rscript "${projectDir}"/bin/differentialExpression.R $gene_expression $params.samplesheet $circRNA_counts $circRNA_raw
        """
    }
}

/*
* DETERMINE miRNA BINDING SITES ON THE PREVIOUSLY DETECTED circRNAs USING miranda
run only if file is not already present
*/
miranda_output = params.out_dir + "/results/binding_sites/output/miranda/bind_sites_raw.out"
if (!file(miranda_output).exists()) {
    miranda_tmp = "${params.out_dir}/results/binding_sites/output/miranda/tmp"
    process miranda {
        label 'process_medium'
        publishDir miranda_tmp, mode: params.publish_dir_mode
        
        input:
        file(circRNA_fasta) from circRNAs_fasta1.splitFasta( by: params.splitter, file: true )

        output:
        file("bind_sites_raw.out") into bind_sites_split

        script:
        """
        miranda $params.mature_fasta $circRNA_fasta -out "bind_sites_raw.out" -quiet
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
    publishDir "${params.out_dir}/results/binding_sites/output/miranda", mode: params.publish_dir_mode
    
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
    publishDir "${params.out_dir}/results/binding_sites/output/miranda", mode: params.publish_dir_mode
    
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
    tarpmir_tmp = "${params.out_dir}/results/binding_sites/output/tarpmir/tmp"
    // RUN TARPMIR ON CHUNKED MRNA FASTAS
    process tarpmir {
        label 'process_medium'
        publishDir tarpmir_tmp, mode: params.publish_dir_mode

        input:
        file(model) from model
        file(mRNA_fasta) from circRNAs_fasta2.splitFasta( by: params.splitter, file: true )

        output:
        file("bindings.bp") into bp_files

        script:
        """
        python3 "${projectDir}/bin/TarPmiR_threading.py" \\
        -a $params.mature_fasta \\
        -b $mRNA_fasta \\
        -m $model \\
        -p $params.p \\
        -t $params.threads \\
        -o "bindings.bp"
        """
    }

    // combine files to one
    bp_files.collectFile(name: "${params.out_dir}/results/binding_sites/output/tarpmir/tarpmir_bp.tsv", newLine: true).into{ tarpmir_bp_file }
    // delete tmp files
    file(tarpmir_tmp).deleteDir()
} else {
    Channel.of( 'null' ).into{ tarpmir_bp_file }
}

/*
 * RUN PITA ANALYSIS FOR circRNAs
 */
if (params.pita) {
    pita_tmp = "${params.out_dir}/results/binding_sites/output/PITA/tmp"
    process PITA {
        label 'process_long'
        publishDir pita_tmp, mode: params.publish_dir_mode
        errorStrategy 'retry'

        input:
        file(circ_fasta) from circRNAs_fasta3.splitFasta( by: params.splitter, file: true )

        output:
        file("circRNA_pita_results.tab") into pita_splits

        script:
        """
        perl $params.pita_path/pita_prediction.pl -utr $circ_fasta -mir $params.mature_fasta -prefix circRNA -l 7-8 -gu 7:0,8:1 -m 7:0,8:0
        """
    }

    // collect all PITA splits
    pita_splits.collectFile(name: "${params.out_dir}/results/binding_sites/output/PITA/circRNA_pita_results.tsv", newLine: true).into{ pita_results }
    // delete tmp directory
    file(pita_tmp).deleteDir()
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
            publishDir "${params.out_dir}/bowtie_index/", mode: params.publish_dir_mode

            input:
            file(fasta) from ch_fasta
            
            output:
            file("${fasta.baseName}*") into ch_generated_bowtie_index
                            
            when: (params.bowtie_index == null)

            script:
            """
            echo "bowtie index is in ${fasta.baseName}"
            bowtie-build $fasta ${fasta.baseName}
            """
        }

        ch_bowtie_index = params.bowtie_index ? Channel.value(file(params.bowtie_index)) : ch_generated_bowtie_index

        /*
        * PERFORM miRNA READ MAPPING USING miRDeep2
        */
        process miRDeep2_mapping {
            label 'process_high'
            publishDir "${params.out_dir}/samples/${sampleID}/miRNA_detection/", mode: params.publish_dir_mode

            input:
            tuple val(sampleID), file(read_file) from ch_smallRNA_reads
            file(index) from ch_bowtie_index.collect()
            file(fasta) from ch_fasta

            output: 
            tuple val(sampleID), file("reads_collapsed.fa"), file("reads_vs_ref.arf") into ch_miRNA_mapping_output

            script:
            """
            gunzip < $read_file > "${sampleID}.fastq"
            mapper.pl "${sampleID}.fastq" -e -h -i -j -k $params.miRNA_adapter -l 18 -m -p ${fasta.baseName} -s "reads_collapsed.fa" -t "reads_vs_ref.arf" -v
            """
        }

        /*
        * PERFORM miRNA QUANTIFICATION USING miRDeep2
        */
        process miRDeep2_quantification {
            label 'process_high'
            publishDir "${params.out_dir}/samples/${sampleID}/miRNA_detection/", mode: params.publish_dir_mode
    
            input:
            tuple val(sampleID), file(reads_collapsed_fa), file(reads_vs_ref_arf) from ch_miRNA_mapping_output
            file(fasta) from ch_fasta

            output:
            file("miRNAs_expressed*") into ch_miRNA_expression_files

            script:
            """
            miRDeep2.pl $reads_collapsed_fa $fasta $reads_vs_ref_arf $params.mature_fasta $params.mature_other_fasta $params.hairpin_fasta -t $params.species -d -v 
            """
        }

        /*
        * MERGE RAW miRNA RESULTS INTO ONE TABLE SUMMARIZING ALL SAMPLES
        */
        process summarize_detected_miRNAs{
            label 'process_medium'

            publishDir "${params.out_dir}/results/miRNA/", mode: params.publish_dir_mode
        
            input:
            file(miRNAs_expressed) from ch_miRNA_expression_files.collect()

            output:
            file("miRNA_counts_raw.tsv") into ch_miRNA_counts_raw

            script:
            """
            Rscript "${projectDir}"/bin/miRNA_summarize_results.R $params.samplesheet $params.out_dir
            """
        }
    }

    /*
    * NORMALIZE RAW miRNA COUNT USING LIBRARY SIZE ESTIMATION
    */
    process normalize_miRNAs{
        label 'process_low'

        publishDir "${params.out_dir}/results/miRNA/", mode: params.publish_dir_mode
        
        input:
        file(miRNA_counts_raw) from ch_miRNA_counts_raw

        output:
        file("miRNA_counts_normalized.tsv") into (ch_miRNA_counts_norm1, ch_miRNA_counts_norm2)

        script:
        """
        Rscript "${projectDir}"/bin/miRNA_results_LibrarySizeEstimation.R $miRNA_counts_raw $params.out_dir
        """
    }

    /*
    * FILTER miRNAs TO REDUCE LOW EXPRESSED ONES
    */
    process filter_miRNAs{
        label 'process_medium'

        publishDir "${params.out_dir}/results/miRNA/", mode: params.publish_dir_mode
        
        input:
        file(miRNA_counts_norm) from ch_miRNA_counts_norm1

        output:
        file("miRNA_counts_filtered.tsv") into (ch_miRNA_counts_filtered1, ch_miRNA_counts_filtered2, ch_miRNA_counts_filtered3)

        script:
        """
        Rscript "${projectDir}"/bin/miRNA_filtering.R $miRNA_counts_norm $params.out_dir $params.sample_percentage $params.read_threshold
        """
    }


    /*
    * FOR ALL POSSIBLE circRNA-miRNA PAIRS COMPUTE PEARSON CORRELATION
    */
    process compute_correlations{
        label 'process_high'

        publishDir "${params.out_dir}/results/sponging/", mode: params.publish_dir_mode
        
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

        publishDir "${params.out_dir}/results/sponging/", mode: params.publish_dir_mode
        
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
        mkdir -p "${params.out_dir}/results/sponging/plots/"
        Rscript "${projectDir}"/bin/correlation_analysis.R $params.samplesheet $miRNA_counts_filtered $circRNA_counts_filtered $correlations $params.out_dir $params.sample_group $miRNA_counts_norm $circRNA_counts_norm
        """
    }

    /*
    * SPONGE ANALYSIS (https://github.com/biomedbigdata/SPONGE)
    */
    if (params.sponge) {
        // USE GIVEN TARGET SYMBOLS OR DEFAULT LOCATED IN DATA
        target_scan_symbols = params.target_scan_symbols ? Channel.value(file(params.target_scan_symbols)) : Channel.value(file(projectDir + "/data/miRNA_target_symbols/hsa_mirWalk_lncbase_21_ENSG.tsv.gz"))
        process SPONGE {
            label 'process_high'
            errorStrategy 'ignore'

            publishDir "${params.out_dir}/results/sponging/SPONGE", mode: params.publish_dir_mode

            input:
            file(gene_expression) from gene_expression
            file(circRNA_counts_filtered) from ch_circRNA_counts_filtered5
            file(mirna_expression) from ch_miRNA_counts_filtered3
            file(miranda_bind_sites) from ch_bindsites_filtered2
            file(target_scan_symbols) from target_scan_symbols
            val(tarpmir) from tarpmir_bp_file
            val(pita) from pita_results

            output:
            file("sponge.RData") into sponge_rimage
            file("plots/*.png") into plots
            file("circRNA/*") into circResults
            file("total/*") into totalResults

            script:
            """
            Rscript "${projectDir}"/bin/SPONGE.R \\
            --gene_expr $gene_expression \\
            --circ_rna $circRNA_counts_filtered \\
            --mirna_expr $mirna_expression \\
            --fdr $params.fdr \\
            --target_scan_symbols $target_scan_symbols \\
            --miranda_data $miranda_bind_sites \\
            --tarpmir_data $tarpmir \\
            --pita_data $pita \\
            --majority_matcher $params.majority_matcher \\
            --normalize
            """
        }
        /*
        * PERFORM SPONGE circRNA DE ANALSIS
        */
        if (params.differential_expression) {
            process SPONGE_DE_ANALYSIS {
                label 'process_low'
                errorStrategy 'ignore'

                publishDir "${params.out_dir}/results/sponging/SPONGE/DE_analysis", mode: params.publish_dir_mode

                input:
                file(sponge_data) from sponge_rimage
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
    }
}
