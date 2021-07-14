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
    tuple val(base), file("${base}.trimmed.R1.fastq.gz")
    file "*"


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
    file "*"

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

process FastQC{ 
container "biocontainers/fastqc:v0.11.9_cv8"
beforeScript 'chmod o+rw .'
cpus 1
publishDir "${params.OUTPUT}/FastQC/", mode: 'symlink'
input: 
    tuple val(base), file(R1) 
output: 
    file "*fastqc*"

script:
"""
#!/bin/bash

echo logging 
ls -lah

fastqc *fastq.gz 

"""
}

process NanoPlot{ 
container "staphb/nanoplot:1.33.0"
beforeScript 'chmod o+rw .'
cpus 2
publishDir "${params.OUTPUT}/NanoPlot/", mode: 'symlink'
input: 
    tuple val(base), file(R1) 
output: 
    file "*.html"

script:
"""
#!/bin/bash

echo logging 
ls -lah

NanoPlot \
    --fastq_rich ${R1} \
    --readtype 2D \
    -t ${task.cpus} \
    -p ${base} \
    --title ${base}

"""
}

process Deduplicate { 
container "broadinstitute/picard"
beforeScript 'chmod o+rw .'
cpus 4
publishDir "${params.OUTPUT}/Deduplicated/", mode: 'symlink'
input: 
    file bam
output: 
    file "*.deduped.bam"
    file "*.txt"

script:
"""
#!/bin/bash

echo logging 
ls -lah

newbase=`echo ${bam} | cut -f1 -d .`

java -jar /usr/picard/picard.jar MarkDuplicates \
    I=${bam} \
    O=\$newbase.deduped.bam \
    M=metrics.txt \
    REMOVE_DUPLICATES=true 
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
    file "${base}.minimap2.sam"

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
    file "${base}.star*"
    file "${base}.starAligned.sortedByCoord.out.bam"
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
publishDir "${params.OUTPUT}/Star_SE/${base}", mode: 'symlink', overwrite: true
container "quay.io/biocontainers/star:2.7.9a--h9ee0642_0"
beforeScript 'chmod o+rw .'
cpus 4
input: 
    tuple val(base), file(r1)
    file starindex
output: 
    file "${base}.star*"
    file "${base}.starAligned.sortedByCoord.out.bam"

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


process Rename { 
//conda "${baseDir}/env/env.yml"
container "mgibio/samtools:1.9"
beforeScript 'chmod o+rw .'
cpus 4
input: 
    file holdFile
output: 
    file "*.sorted.bam"
script:
"""
#!/bin/bash

#logging
echo "ls of directory" 
ls -lah 

count=`ls -1 *.sam 2>/dev/null | wc -l`
if [ \$count != 0 ]
then 
echo "running samtools for nanopore file"
for i in *.sam
do
base=`basename -s ".bam" \$i` 
samtools view -Sb -@ ${task.cpus} \$i > \$base.bam
done
fi 

for i in *.bam
do
echo sorting \$i
newbase=`echo \$i | cut -f1 -d .`
samtools sort \$i -@ ${task.cpus} -o \$newbase.sorted.bam
done

"""
}


process Chromoplots { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/generate_chromo_plots/", mode: 'symlink'
container "quay.io/vpeddu/exote:latest"
beforeScript 'chmod o+rw .'
cpus 2
errorStrategy 'retry'
maxRetries 3

input: 
    file bam
    file repeat_GTF
    file genomic_GTF
    file annotations
    file chr_info
    val NANOPORE
    val PAIRED_END

output: 
    tuple file("*.csv"), file('*.rData'), file("*.html")

script:
"""
#!/bin/bash

echo "conda activate TE_Flow" >> ~/.bashrc
source ~/.bashrc

#logging
echo "ls of directory" 
ls -lah 

echo "unzipping"
gunzip *.gz
ls -lah 


echo paired ${PAIRED_END}
echo nanopore ${NANOPORE}

Rscript --vanilla ${baseDir}/bin/chromo_plots.r \
    ${chr_info} \
    ${bam} \
    ${repeat_GTF} \
    ${annotations} \
    ${task.cpus} \
    `basename -s ".bam" *.bam` \
    ${genomic_GTF} \
    ${NANOPORE} \
    ${PAIRED_END}


"""
}

//could parallelize later but reasonably fast now
process Biotyping{ 
container "quay.io/vpeddu/alfa"
beforeScript 'chmod  o+rw .'
cpus 8
publishDir "${params.OUTPUT}/Alfa/", mode: 'symlink'
input: 
    file bam
    file genomic_GTF
output: 
    file "*.pdf"

shell:
'''
#!/bin/bash
echo logging ls
ls -lah


#sorting GTF
sort -k1,1 -k4,4n -k5,5nr !{genomic_GTF} > sorted_gtf.gtf
echo "finished sort"

#renaming biotype field 
sed -i 's/gene_type/gene_biotype/g' sorted_gtf.gtf
echo "finished rename"

#index generation
alfa -a sorted_gtf.gtf -g index.alfaindex -p !{task.cpus}
echo "finished index generation"

for i in *.bam
do
    base=`basename $i ".bam"` 
    echo working on $base
    alfa -g index.alfaindex \
    --bam $i $base \
    -p !{task.cpus} \
    -d 3
    mv ALFA_plots.Biotypes.pdf $base.Biotypes.pdf
    mv ALFA_plots.Categories.pdf $base.Categories.pdf
done 

'''
}

process TEtranscripts{ 
container "quay.io/vpeddu/alfa"
beforeScript 'chmod o+rw .'
cpus 1
publishDir "${params.OUTPUT}/TEtranscripts/", mode: 'symlink'
input: 
    file bam
    file genomic_GTF
    file repeat_GTF
output: 
    file "*.TEtranscripts"

script:
"""
#!/bin/bash

echo logging 
ls -lah

base=`basename \$i` .bam

echo base is \$base

TEcount \
    --format BAM \
    --mode uniq \
    --sortByPos \
    -b ${bam} \
    --GTF ${genomic_GTF} \
    --TE ${repeat_GTF}

mkdir \$base.TEtranscripts
mv * \$base.TEtranscripts
"""
}

process MultiQC{ 
container "ewels/multiqc:1.10.1"
beforeScript 'chmod o+rw .'
cpus 1
publishDir "${params.OUTPUT}/MultiQC/", mode: 'symlink'
input: 
    file infolder
    file renamehold
output: 
    file "*"

script:
"""
#!/bin/bash

echo logging 
ls -lah

ls /usr/local/bin/
multiqc .

echo mutliqc_finish > done.txt

"""
}