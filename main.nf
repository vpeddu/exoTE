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
include { Minimap2 }from './modules.nf'
include { Star_PE }from './modules.nf'
include { Star_SE }from './modules.nf'
include { Rename }from './modules.nf'
include { Chromoplots }from './modules.nf'

// these params must be defined in nextflow.config? or is it defined on the fly via the CLI?
if (!params.INPUT_FOLDER){ exit 1, "Must provide folder containing input files with --CONTROL_INPUT" }
//if (!params.ADAPTER_CHOICE){ exit 1, "Must provide adapter choice with --ADAPTER_CHOICE (intra, exo)" }
//if (!params.SALMON_INDEX){ exit 1, "Must provide adapter choice with --SALMON_INDEX" }
//if (!params.OUTPUT){ exit 1, "Must provide adapter choice with --OUTPUT" }
//if (!params.FASTA){ exit 1, "Must provide reference fasta with --FASTA" }
//if (!params.GTF){ exit 1, "Must provide reference GTF with --GTF" }



	//.ifEmpty { exit 1, "Star index not found: ${params.STAR_INDEX}" }

repeat_GTF= file(params.REPEAT_GTF)
genomic_GTF=file(params.GENOMIC_GTF)

workflow{
    if ( params.NANOPORE ) { 
        Minimap2_ref = file( params.MINIMAP2_REFERENCE )
        input_read_Ch = Channel
            .fromPath("${params.INPUT_FOLDER}*.fastq.gz")
            .map { it -> [it.name.replace(".fastq.gz", ""), file(it)]}
        Minimap2( 
            input_read_Ch,
            Minimap2_ref
        )
        }
    // defined at CLI
    else if ( params.ILLUMINA){ 
        Star_index_Ch = Channel
            .fromPath(params.STAR_INDEX)
            
        if ( params.PAIRED_END ){
        input_read_Ch = Channel
            .fromFilePairs("${params.INPUT_FOLDER}**_R{1,2}*.fastq.gz")
            .map { it -> [it[0], it[1][0], it[1][1]]}

        Fastp_PE(
            input_read_Ch
            )
        Star_PE(
            // access output of preceeding process
            Fastp_PE.out,
            // collects all items emitted by a channel to a list, return
            Star_index_Ch.collect()
            )
        }

        // start single end workflow
        else if ( params.SINGLE_END ) { 
            input_read_Ch = Channel
                .fromPath("${params.INPUT_FOLDER}*.fastq.gz")
                .map { it -> [it.name.replace(".fastq.gz", ""), file(it)]}

            Fastp_SE(
                input_read_Ch
                )
            Star_SE(
                Fastp_SE.out,
                Star_index_Ch.collect()
                )
            }
        }

// combining paired and single end workflows here
// converts minimap2 sam to bam here too
    if (params.NANOPORE){ 
        Rename( 
            Minimap2.out.collect()
        )
    }
    else if ( params.PAIRED_END ){ 
        Rename(
        Star_PE.out.toList()
        )
    } 
    else if ( params.SINGLE_END) { 
        Rename(
        Star_SE.out.collect()
        )
    }
    Chromoplots( 
        Rename.out,
        file(repeat_GTF),
        file(genomic_GTF)

    )
}