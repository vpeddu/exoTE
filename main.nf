#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def helpMessage() {
// log.info "foo" prints to STDOUT and nextflow.log
    log.info"""
    TE_Flow: Pipeline that combines standard RNA-seq quantification approaches with TE-aware
    steps to generate a comprehensive accounting of transcriptomic activity.
    Usage:
    nextflow run rreggiar/TE_Flow \\
        --INPUT_FOLDER test_fastq_paired_small/ 
        --ADAPTER_CHOICE intra 
        --PAIRED_END True 
        --OUTPUT testout 
        --SALMON_INDEX test/test_salmon_index 
        -with-docker Ubuntu:18.04
    """.stripIndent()
}

// show help message
// accessible via nextflow run $project --help
params.help = false
// The params scope allows you to define parameters that will be accessible in the pipeline script
// my guess is params.help is a global boolean that depends on --help at CLI
if (params.help){
    // if help == T
    helpMessage()
    // run helpMessage()
    exit 0
    // clean exit
}


include { Fastp_SE }from './modules.nf'
include { Fastp_PE }from './modules.nf'
include { Salmon_PE }from './modules.nf'
include { Salmon_SE }from './modules.nf'
include { Salmon_hold }from './modules.nf'
include { Generate_col_data }from './modules.nf'

// these params must be defined in nextflow.config? or is it defined on the fly via the CLI?
if (!params.INPUT_FOLDER){ exit 1, "Must provide folder containing input files with --CONTROL_INPUT" }
//if (!params.ADAPTER_CHOICE){ exit 1, "Must provide adapter choice with --ADAPTER_CHOICE (intra, exo)" }
//if (!params.SALMON_INDEX){ exit 1, "Must provide adapter choice with --SALMON_INDEX" }
if (!params.OUTPUT){ exit 1, "Must provide adapter choice with --OUTPUT" }
if (!params.FASTA){ exit 1, "Must provide reference fasta with --FASTA" }
if (!params.GTF){ exit 1, "Must provide reference GTF with --GTF" }

// adapters = [
//     file("${baseDir}/adapters/NexteraPE-PE.fa"),
//     file("${baseDir}/adapters/TruSeq3-PE.fa")
]

// salmon_index_Ch = Channel
// 	.fromPath(params.SALMON_INDEX)
// 	.ifEmpty { exit 1, "Salmon index not found: ${params.SALMON_INDEX}" }



workflow{
    // defined at CLI
    if ( params.PAIRED_END ){
        input_read_Ch = Channel
            .fromFilePairs("${params.INPUT_FOLDER}**_R{1,2}*.fastq.gz")
            // little confused about what's going on here, mapping from ^
            // to ... a list? what is `it`?
            .map { it -> [it[0], it[1][0], it[1][1]]}

        Trimmomatic_PE(
            input_read_Ch,
            params.ADAPTER_CHOICE, 
            adapters
            )
        Salmon_PE(
            // access output of preceeding process
            Trimmomatic_PE.out,
            // collects all items emitted by a channel to a list, return
            salmon_index_Ch.collect()
            )
        }

    // start single end workflow
    else{ 
        input_read_Ch = Channel
            .fromPath("${params.INPUT_FOLDER}*.fastq.gz")
            .map { it -> [it.name.replace(".fastq.gz", ""), file(it)]}

        Trimmomatic_SE(
            input_read_Ch,
            params.ADAPTER_CHOICE, 
            adapters
            )
        Salmon_SE(
            Trimmomatic_SE.out,
            salmon_index_Ch.collect()
            )
        }
    // combining paired and single end workflows here
    if ( params.PAIRED_END){ 
        Salmon_hold(
        Salmon_PE.out.toList()()
        )
    } else { 
        Salmon_hold(
        Salmon_SE.out.collect()
        )
    }
    Generate_col_data( 
        file(params.FASTA),
        file(params.GTF),
        Salmon_hold.out,
        salmon_index_Ch.collect(),
        file(params.TXTG)

    )
}