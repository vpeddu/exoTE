library(chromoMap)
library(Rsamtools)
library(Rsubread)
library(tidyverse)
library(biomaRt)
library(htmlwidgets)
library(readr)

print('libraries imported')


args = commandArgs(trailingOnly=TRUE)

chr_info<-read.csv(args[1], header = FALSE)
bam<-as.data.frame(scanBam(args[2]))
print('read bam')
repeat_annotations<-as.data.frame(rtracklayer::import(args[3]))
print('read annotations')

base= args[6]

# file was generated on 6/16/21 from table browser rmsk hg38 full annotations
# subsetted to only include LINE and SINE repclass
# each line is chromosome.start.end.strand.element_name
full_annotations<-read_csv(args[3])

repeat_annotations<-repeat_annotations[repeat_annotations$gene_id %in% full_annotations,]

annotations<-repeat_annotations[,c(10,1,2,3)]

# data cleaning 
annotations<-annotations[-which(grepl('alt',annotations$seqnames)),]
annotations<-annotations[-which(grepl('random',annotations$seqnames)),]
annotations<-annotations[-which(grepl('chrUn',annotations$seqnames)),]
annotations$gene_id<-as.character(annotations$gene_id)

#annotations<-head(annotations,50000)

#colnames(annotations)[5]<-'strand'

chr_info<-chr_info[which(chr_info$V1 %in% annotations$seqnames),]
chr_info<-chr_info[order(chr_info$V1),]
chr_info<-chr_info[,c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr11","chr10","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")]
#chr_info<-head(chr_info)

print('files read in. starting featurecounts ')

# featurecounts with rmsk gtf to assign reads 
TE_fc_out<-featureCounts(files = args[2],
                         annot.ext = args[3],
                         isGTFAnnotationFile = T,
                         nthreads = args[5])

print('TE featurecounts done ')
# "genomic" genes
genomic_fc_out<-featureCounts(files = args[2],
                              annot.ext = args[4],
                              isGTFAnnotationFile = T,
                              nthreads = args[5])

print('Genomic featurecounts done')

TE_fc_counts<-as.data.frame(TE_fc_out$counts)
TE_fc_counts<-rownames_to_column(TE_fc_counts, var = 'repeat_element')
TE_fc_counts<-TE_fc_counts[which(TE_fc_counts$repeat_element %in% annotations$gene_id),]
colnames(TE_fc_counts)[2]<-'counts'
TE_fc_counts<-TE_fc_counts[TE_fc_counts$counts > 2,]


TE_chromo_df <- data.frame(TE_fc_counts$repeat_element,
                           sapply(strsplit(as.character(TE_fc_counts$repeat_element), "[.]"), `[`, 1),
                           sapply(strsplit(as.character(TE_fc_counts$repeat_element), "[.]"), `[`, 2),
                           sapply(strsplit(as.character(TE_fc_counts$repeat_element), "[.]"), `[`, 3),
                           TE_fc_counts$counts
)


genomic_fc_counts<-as.data.frame(genomic_fc_out$counts)


print('starting biomaRt')
# filling in gene annotations for genomic with biomart
mart <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart)
gene_lookup<-getBM(attributes = c("ensembl_gene_id_version",'start_position','end_position','chromosome_name'),
                   mart = ensembl)
gene_lookup$chromosome_name<-paste0('chr',gene_lookup$chromosome_name)
gene_lookup$chromosome_name[gene_lookup$chromosome_name == 'chrMT'] <- 'chrM'
colnames(gene_lookup)[1]<-'gene'

genomic_fc_counts<-rownames_to_column(genomic_fc_counts, var = 'gene')
colnames(genomic_fc_counts)[2]<-'counts'
genomic_fc_counts<-genomic_fc_counts[genomic_fc_counts$counts > 2,]

genomic_chromo_df <- dplyr::full_join(genomic_fc_counts,gene_lookup, by = 'gene')
genomic_chromo_df <- genomic_chromo_df[,c(1,5,3,4,2)]
### pileup test
# pparam<-PileupParam(min_nucleotide_depth=3, 
#             distinguish_strands=TRUE)

# pup<-pileup('/Users/vikas/Documents/UCSC/lab/kim_lab/exosomal_rna/ctrl.10.0.4.bam', 
#             index='/Users/vikas/Documents/UCSC/lab/kim_lab/exosomal_rna/ctrl.10.0.4.bam.bai',
#             pileupParam = pparam)
#test<-chromo_df[1:100,]
#chromo_df<-chromo_df[chromo_df$fc_counts.repeat_element > 0,]
#test_chr<-head(chr_info,1)
#colnames(test)<-c('a','b','c','d')
#test$a<-c('1','2','3','4','5','6')
print("saving image")
save.image(paste0(args[6],'chromo_plots.rData'))


print('making plots')
repeat_plot<-chromoMap(list(chr_info),list(TE_chromo_df),
          data_based_color_map = T,
          plot_height = 10,
          title = paste(args[6],'repeats'), 
          canvas_height = 1000000,
          data_type = "numeric",
          #plots = "bar",
          heat_map = F,
          ch_gap = .1,
          # canvas_width = 600,
          top_margin = 1)
# left_margin = 25)


genomic_plot<-chromoMap(list(chr_info),list(genomic_chromo_df),
          data_based_color_map = T,
          plot_height = 10,
          canvas_height = 1000000,
          title = paste(args[6],'genomic'), 
          data_type = "numeric",
          #plots = "bar",
          heat_map = T,
          ch_gap = .1,
          # canvas_width = 600,
          top_margin = 1)

saveWidget(
  repeat_plot,
  paste0(args[6],'_repeat_chromoplot.html')
)
saveWidget(
  genomic_plot,
  paste0(args[6],'_genomic_chromoplot.html')
)

write.csv(genomic_chromo_df,'genomic_chromo_df.csv')
write.csv(TE_chromo_df,'TE_chromo_df.csv')
