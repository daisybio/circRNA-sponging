include { STAR_ALIGN } from '../modules/nf-core/star/align/main.nf'
include { CIRCEXPLORER2_PARSE } from '../modules/nf-core/circexplorer2/parse/main'
include { CIRCEXPLORER2_ANNOTATE } from '../modules/nf-core/circexplorer2/annotate/main'
include { SUMMARIZE } from '../modules/local/summarize.nf'
include { BEDTOOLS_GETFASTA } from '../modules/nf-core/bedtools/getfasta/main.nf'

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
        array = [ [id: row.sample, single_end: params.single_end], [ file(row.totalRNA1), file(row.totalRNA2) ] ]
    }
    if (params.differential_expression){
        if (!row.containsKey("condition")) {
            exit 1, "Error: Condition information missing!"
        }
    }
    return array
}

ch_totalRNA_reads = Channel.fromPath(params.samplesheet) // [meta, [reads1, reads2]]
   .splitCsv( header:true, sep:'\t')
   .map { get_circRNA_paths(it) }

workflow CIRCRNA {
    take:
        ch_star_index
        ch_gtf
        ch_genome_fasta
        ch_bed12
    
    main:
        STAR_ALIGN (
            ch_totalRNA_reads,
            ch_star_index.map{ [[], it] },
            ch_gtf.map{ [[], it] },
            true,
            false,
            false
        )

        CIRCEXPLORER2_PARSE (
            STAR_ALIGN.out.junction
        )

        CIRCEXPLORER2_ANNOTATE (
            CIRCEXPLORER2_PARSE.out.junction,
            ch_genome_fasta,
            ch_bed12
        )

        SUMMARIZE (CIRCEXPLORER2_ANNOTATE.out.txt.map{ it[1] }.collect())

        BEDTOOLS_GETFASTA (
            SUMMARIZE.out,
            ch_genome_fasta
        )
}