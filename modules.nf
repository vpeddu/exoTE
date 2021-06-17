// a process is the basic processing `primitive` to exec a user script

// nextflow process syntax:

// process < name > {
//  [ directives ]
//  
//  input:
//    < process inputs >
//
//  output:
//    < process outputs >
//
//  when:
//    < condition >
//
//  [script|shell|exec]:
//  < user script to be executed >
//}


process Fastp_PE { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/fastp_PE/${base}", mode: 'symlink', overwrite: true
container "bromberglab/fastp"
beforeScript 'chmod o+rw .'
cpus 4
input: 
    tuple val(base), file(r1), file(r2)
    val adapter
    tuple file(a1), file(a2)
output: 
    tuple val(base), file("${base}.trimmed.R1.fastq.gz"), file("${base}.trimmed.R2.fastq.gz")

script:
"""
#!/bin/bash

#logging
echo "ls of directory" 
ls -lah 

echo "running fastp on ${base}"

fastp -w ${task.cpus} \
    -i ${r1} \
    -I ${r2} \
    -o ${base}.trimmed.R1.fastq.gz \
    -O ${base}.trimmed.R2.fastq.gz
"""
}

process Fastp_SE { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/Trimmomatic_SE/${base}", mode: 'symlink'
container "staphb/trimmomatic"
beforeScript 'chmod o+rw .'
cpus 4
input: 
    tuple val(base), file(r1)
    val adapter
    tuple file(a1), file(a2)
output: 
    tuple val(base), file("${base}.trimmed.R1.fastq.gz")

script:
"""
#!/bin/bash

#logging
echo "ls of directory" 
ls -lah 

echo "running fastp on ${base}"

fastp \
    -w ${task.cpus} \
    -i ${r1} \
    -o ${base}.trimmed.R1.fastq.gz 
"""
}

process Salmon_PE { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/Salmon_PE/${base}", mode: 'symlink', overwrite: true
container "combinelab/salmon"
beforeScript 'chmod o+rw .'
cpus 4
input: 
    tuple val(base), file(r1), file(r2)
    file salmonIndex
output: 
    tuple val(base), file("quant.sf")
script:
"""
#!/bin/bash

#logging
echo "ls of directory" 
ls -lah 

echo "salmon version:" \$(salmon --version)

    salmon quant \
        -i ${salmonIndex} \
        --libType A \
        -1 ${r1} \
        -2 ${r2} \
        -p ${task.cpus} \
        --validateMappings \
        --gcBias \
        --seqBias \
        --recoverOrphans \
        --rangeFactorizationBins 4 \
        --output . 

"""
}
process Salmon_SE { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/Salmon_SE/", mode: 'symlink'
container "combinelab/salmon"
beforeScript 'chmod o+rw .'
cpus 4
input: 
    tuple val(base), file(r1)
    file salmonIndex
output: 
    file "${base}"
script:
"""
#!/bin/bash

#logging
echo "ls of directory" 
ls -lah 

echo "salmon version:" \$(salmon --version)

    salmon quant \
        -i ${salmonIndex} \
        --libType A \
        -r ${r1} \
        -p ${task.cpus} \
        --validateMappings \
        --gcBias \
        --seqBias \
        --recoverOrphans \
        --rangeFactorizationBins 4 \
        --output . 

mkdir ${base}
mv quant.sf ${base}
"""
}


process Salmon_hold { 
//conda "${baseDir}/env/env.yml"
container "combinelab/salmon"
beforeScript 'chmod o+rw .'
cpus 1
input: 
    file quantfolder
output: 
    file "salmon_files"
script:
"""
#!/bin/bash

#logging
echo "ls of directory" 
ls -lah 

mkdir salmon_files
ls | grep -v "salmon_files" | xargs -I % mv % salmon_files
"""
}


process Generate_col_data { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/generate_col_data/", mode: 'symlink'
//container "quay.io/vpeddu/te_flow"
container "testing"
beforeScript 'chmod o+rw .'
cpus 4

input: 
 file FASTA
 file GTF
 file SALMON_FILES
 file SALMON_INDEX
 file TXTG

output: 
    tuple file("results/assays.h5"), file('results/se.rds')

script:
"""
#!/bin/bash

echo "conda activate TE_Flow" >> ~/.bashrc
source ~/.bashrc

#logging
echo "ls of directory" 
ls -lah 

# run makeLinkedTxome 
Rscript --vanilla ${baseDir}/R/generateLinkedTxome.R ${SALMON_INDEX}/ ${FASTA} ${GTF}

echo "built json "

Rscript --vanilla ${baseDir}/R/generateColData.R . txome.json shit lol output/ ${TXTG}


"""
}