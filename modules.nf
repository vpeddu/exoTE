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
publishDir "${params.OUTPUT}/Fastp_SE/${base}", mode: 'symlink'
container "bromberglab/fastp"
beforeScript 'chmod o+rw .'
cpus 4
input: 
    tuple val(base), file(r1)
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

process Minimap2 { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/Minimap2/${base}", mode: 'symlink'
container "biocontainers/minimap2:v2.15dfsg-1-deb_cv1"
beforeScript 'chmod o+rw .'
cpus 4
input: 
    tuple val(base), file(r1)
    file(Minimap2_ref)
output: 
    tuple val(base), file("${base}.minimap2.sam")

script:
"""
#!/bin/bash

#logging
echo "ls of directory" 
ls -lah 

echo "running Minimap2 on ${base}"

minimap2 \
    -ax map-ont \
    -t ${task.cpus} \
    ${Minimap2_ref} \
    ${r1} > \
    ${base}.minimap2.sam
"""
}

process Star_PE { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/Star_PE/${base}", mode: 'symlink', overwrite: true
container "quay.io/biocontainers/star:2.7.9a--h9ee0642_0"
beforeScript 'chmod o+rw .'
cpus 4
input: 
    tuple val(base), file(r1), file(r2)
    file starindex
output: 
    tuple val(base), file("${base}.starAligned.sortedByCoord.out.bam")
script:
"""
#!/bin/bash

#logging
echo "ls of directory" 
ls -lah 

STAR   \
    --runThreadN ${task.cpus}  \
    --genomeDir ${starindex}   \
    --readFilesIn ${r1} ${r2} \
    --readFilesCommand zcat      \
    --outFileNamePrefix ${base}.star   \
    --outSAMtype BAM   SortedByCoordinate   
"""
}
process Star_SE { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/STAR_SE/", mode: 'symlink'
container "quay.io/biocontainers/star:2.7.9a--h9ee0642_0"
beforeScript 'chmod o+rw .'
cpus 4
input: 
    tuple val(base), file(r1)
    file starindex
output: 
    tuple val(base), file("${base}.starAligned.sortedByCoord.out.bam")
script:
"""
#!/bin/bash

#logging
echo "ls of directory" 
ls -lah 

STAR   \
    --runThreadN ${task.cpus}  \
    --genomeDir ${starindex}   \
    --readFilesIn ${r1}  \
    --readFilesCommand zcat      \
    --outFileNamePrefix ${base}.star   \
    --outSAMtype BAM   SortedByCoordinate   
"""
}


process Hold { 
//conda "${baseDir}/env/env.yml"
container "mgibio/samtools:1.9"
beforeScript 'chmod o+rw .'
cpus 1
input: 
    tuple val(base), file(holdFile)
output: 
    file "*.bam"
script:
"""
#!/bin/bash

#logging
echo "ls of directory" 
ls -lah 

count=`ls -1 *.sam 2>/dev/null | wc -l`
if [ $count != 0 ]
then 
echo "running samtools for nanopore file"
for i in *.sam
do
samtools view -Sb -@ ${task.cpus} $i >  
fi 


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