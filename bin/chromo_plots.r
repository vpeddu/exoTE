library(chromoMap)
library(Rsamtools)
library(Rsubread)
library(tidyverse)
library(biomaRt)
library(htmlwidgets)
library(readr)
library(karyoploteR)
#library(chromstaR)

print('libraries imported')

args = commandArgs(trailingOnly=TRUE)

chr_info<-read.csv(args[1], header = FALSE)
bam<-as.data.frame(scanBam(args[2]))
print('read bam')
repeat_annotations<-as.data.frame(rtracklayer::import(args[3]))
print('read annotations')

base= args[6]

is_nanopore_reads = if (args[8] == 'true') TRUE else FALSE
is_paired_end = if (args[9] == 'true') TRUE else FALSE


# file was generated on 6/16/21 from table browser rmsk hg38 full annotations
# subsetted to only include LINE and SINE repclass
# each line is chromosome.start.end.strand.element_name
full_annotations<-read_csv(args[4])
print('read csv')

repeat_annotations<-repeat_annotations[repeat_annotations$gene_id %in% full_annotations,]

annotations<-repeat_annotations[,c(10,1,2,3)]


# data cleaning 
annotations<-annotations[-which(grepl('alt',annotations$seqnames)),]
annotations<-annotations[-which(grepl('random',annotations$seqnames)),]
annotations<-annotations[-which(grepl('chrUn',annotations$seqnames)),]
annotations<-annotations[-which(grepl('_random',annotations$seqnames)),]
annotations$gene_id<-as.character(annotations$gene_id)

#annotations<-head(annotations,50000)

#colnames(annotations)[5]<-'strand'

#chr_info<-chr_info[which(chr_info$V1 %in% annotations$seqnames),]
chr_info<-chr_info[order(chr_info$V1),]
#chr_info<-chr_info[,c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr11","chr10","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")]
#chr_info<-head(chr_info)

print('files read in. starting featurecounts ')

# featurecounts with rmsk gtf to assign reads 
TE_fc_out<-featureCounts(files = args[2],
                         annot.ext = args[3],
                         isGTFAnnotationFile = T,
                         nthreads = args[5],
                         isLongRead = is_nanopore_reads,
                         isPairedEnd = is_paired_end
                         )

print('TE featurecounts done ')
# "genomic" genes
genomic_fc_out<-featureCounts(files = args[2],
                              annot.ext = args[7],
                              isGTFAnnotationFile = T,
                              nthreads = args[5],
                              isLongRead = is_nanopore_reads,
                              isPairedEnd = is_paired_end)

print('Genomic featurecounts done')

TE_fc_counts<-as.data.frame(TE_fc_out$counts)
TE_fc_counts<-rownames_to_column(TE_fc_counts, var = 'repeat_element')
#TE_fc_counts<-TE_fc_counts[which(TE_fc_counts$repeat_element %in% annotations$gene_id),]
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
#genomic_fc_counts<-genomic_fc_counts[order(genomic_fc_counts$counts, decreasing = TRUE),]
#genomic_fc_counts<-genomic_fc_counts[1:1e2,]

genomic_chromo_df <- dplyr::full_join(genomic_fc_counts,gene_lookup, by = 'gene')
genomic_chromo_df <- genomic_chromo_df[,c(1,5,3,4,2)]
genomic_chromo_df<-genomic_chromo_df[complete.cases(genomic_chromo_df),]


TE_chromo_df<-TE_chromo_df[-which(grepl('alt',TE_chromo_df$sapply.strsplit.as.character.TE_fc_counts.repeat_element...........)),]
TE_chromo_df<-TE_chromo_df[-which(grepl('random',TE_chromo_df$sapply.strsplit.as.character.TE_fc_counts.repeat_element...........)),]
TE_chromo_df<-TE_chromo_df[-which(grepl('chrUn',TE_chromo_df$sapply.strsplit.as.character.TE_fc_counts.repeat_element...........)),]

#genomic_chromo_df<-genomic_chromo_df[-which(grepl('alt',genomic_chromo_df$chromosome_name)),]
#genomic_chromo_df<-genomic_chromo_df[-which(grepl('random',genomic_chromo_df$chromosome_name)),]
#genomic_chromo_df<-genomic_chromo_df[-which(grepl('chrUn',genomic_chromo_df$chromosome_name)),]


TE_chromo_df<-TE_chromo_df[-which(grepl('\\(',TE_chromo_df$TE_fc_counts.repeat_element)),]
TE_chromo_df<-TE_chromo_df[rep(row.names(TE_chromo_df),TE_chromo_df$TE_fc_counts.counts),]
genomic_chromo_df<-genomic_chromo_df[rep(row.names(genomic_chromo_df),genomic_chromo_df$counts),]

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
          heat_map = F,
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

write.csv(genomic_chromo_df,paste0(args[6],'_genomic_chromo_df.csv'))
write.csv(TE_chromo_df,paste0(args[6],'_TE_chromo_df.csv'))




###TESTING
# kp <- plotKaryotype(genome = "hg38", chromosomes = "chr4")
# kpAddBaseNumbers(kp, tick.dist = 50000, add.units = TRUE)
# kp <- kpPlotBAMCoverage(kp, data = file.path('/Users/vikas/Documents/UCSC/lab/kim_lab/exosomal_rna/chromo_plots/testing/BioIVT_panc6_2-5_S88.starAligned.sortedByCoord.out.bam'))

# kp <- plotKaryotype(genome = "hg38", chromosomes = "chr1")
# kpAddBaseNumbers(kp, tick.dist = 50000, add.units = TRUE)
# kp <- kpPlotBAMCoverage(kp, 
#                         data=file.path('/Users/vikas/Documents/UCSC/lab/kim_lab/exosomal_rna/chromo_plots/testing/BioIVT_panc6_2-5_S88.starAligned.sortedByCoord.out.bam'), 
#                         max.valid.region.size = 2e8,
#                         col="gold", 
#                         chromosomes = 'chr1',
#                         border="red",
#                         density=20
#                         )

## AHHHHHHHHHHHHHHHHH 
# bed_in<-rtracklayer::import.bed('/Users/vikas/Documents/UCSC/lab/kim_lab/exosomal_rna/chromo_plots/testing/test.bed',
#                                 genome = 'hg38',
#                                 )
# gr.cov <- as(coverage(bed_in), "GRanges")

# autoplot('/Users/vikas/Documents/UCSC/lab/kim_lab/exosomal_rna/chromo_plots/testing/BioIVT_panc6_2-5_S88.starAligned.sortedByCoord.out.bam',
#          method = "estimate", 
#          which = paste0("chr", 1:22), aes(y = log(..coverage..)))

# bam_in<-readGAlignmentPairs('/Users/vikas/Documents/UCSC/lab/kim_lab/exosomal_rna/chromo_plots/testing/BioIVT_panc6_2-5_S88.starAligned.sortedByCoord.out.bam',
#                             use.names = TRUE)
# #dn$pvalue <- runif(length(dn)) * 10
# p <- autoplot(seqinfo(bam_in)) + layout_karyogram(bam_in, aes(x = start, y = width),
#                                               geom = "point", color = "#fdc086")

# chroms<-c('chr1',
#           'chr2',
#           'chr3',
#           'chr4',
#           'chr5',
#           'chr6',
#           'chr7',
#           'chr8',
#           'chr9',
#           'chr10',
#           'chr11',
#           'chr12',
#           'chr13',
#           'chr14',
#           'chr15',
#           'chr16',
#           'chr17',
#           'chr18',
#           'chr19',
#           'chr20',
#           'chr21',
#           'chr22',
#           'chrX',
#           'chrY',
#           'chrM')
# bed_in<-readBamFileAsGRanges('/Users/vikas/Documents/UCSC/lab/kim_lab/exosomal_rna/chromo_plots/testing/BioIVT_panc6_2-5_S88.starAligned.sortedByCoord.out.bam',
#                              chromosomes = chroms,
#                              pairedEndReads = TRUE)

# x<-GRcoverage(Object=bed_in, bam='/Users/vikas/Documents/UCSC/lab/kim_lab/exosomal_rna/chromo_plots/testing/BioIVT_panc6_2-5_S88.starAligned.sortedByCoord.out.bam', Nnorm=F, Snorm=F)


# start.time <- Sys.time()
# kpPlotBAMDensity(kp, data='/Users/vikas/Documents/UCSC/lab/kim_lab/exosomal_rna/chromo_plots/testing/gencode.test.bam',
#                  col="gold", 
#                  border="blue",
#                  density=2000)

# kpPlotBAMDensity(kp, data='/Users/vikas/Documents/UCSC/lab/kim_lab/exosomal_rna/chromo_plots/testing/TE.test.bam',
#                  col="gold", 
#                  border="red",
#                  density=2000)


# data.points.colors <- c("#FFBD07", "#00A6ED",  "#FF3366", "#8EA604", "#C200FB")
# kp <- plotKaryotype(genome = 'hg38',
#                     plot.type = 1,
#                     chromosomes = c("chr14", "chr2", "chr3"))
# kp <- kpPlotBAMDensity(kp, data='/Users/vikas/Documents/UCSC/lab/kim_lab/exosomal_rna/chromo_plots/testing/gencode_in_regions.bam', 
#                        col=data.points.colors[1], 
#                        border=NA, 
#                        r0=0.5, 
#                        r1=0,
#                        #density = 10000,
#                        #window.size = 1e3,
#                        normalize = FALSE)
# kpAxis(kp, r0=0.5, r1=0, ymax=15, numticks = 1)
# #kpAbline(kp, h=0.4, data.panel = 2, r0=0.2, r1=0, col=data.points.colors[3])

# kp <- kpPlotBAMDensity(kp, data='/Users/vikas/Documents/UCSC/lab/kim_lab/exosomal_rna/chromo_plots/testing/TE_in_regions.bam', 
#                        col=data.points.colors[2], 
#                        border=NA, 
#                        r0=0.5, 
#                        r1=1,
#                        #density = 10000,
#                        #window.size = 1e3,
#                        normalize = FALSE,)
# kpAxis(kp, r0=0.5, r1=1, ymax=15, numticks = 1)


# end.time <- Sys.time()





# ####
# kp <- plotKaryotype(genome = 'hg38', plot.type = 2, chromosomes = c("chr1" ))# , "chr2", "chr3"))
# data.points.colors <- c("#FFBD07", "#00A6ED",  "#FF3366", "#8EA604", "#C200FB")


# kpAxis(kp, ymin = 0, ymax = 1, r0=0, r1=0.8, numticks = 5, col="#666666", cex=0.5)
# kpPoints(kp, data=bed_in, pch=16, cex=0.5, col=dp.colors, r0=0, r1=0.8)

# #Mean and sd of the data points.  
# for(chr in seqlevels(kp$genome)) {
#   chr.dp <- sort(keepSeqlevels(x = bed_in, value = chr, pruning.mode = "coarse"))
#   rmean <- rollmean(chr.dp$n, k = 6, align = "center")  
#   rsd <- rollapply(data = chr.dp$y, FUN=sd, width=6)
#   kpLines(kp, chr = chr, x=start(chr.dp)[3:(length(chr.dp)-3)], y=rmean, col=data.points.colors[3], r0=0, r1=0.8)
#   kpPlotRibbon(kp, chr=chr, data=chr.dp[3:(length(chr.dp)-3)], y0=rmean-rsd, y1=rmean+rsd, r0=0, r1=0.8, col="#FF336633", border=NA)
# }
# ### Data Panel 2 ###

# #medium regions and their coverage

# kpPlotRegions(kp, data = bed_in, r0 = 0.2, r1=1, border=NA, data.panel=2)
# kpPlotCoverage(kp, data=bed_in, r0=0.2, r1=0, col=data.points.colors[2], data.panel = 2)
# kpPlotCoverage(kp, data=bed_in, r0=0.2, r1=0.12, col=data.points.colors[1], data.panel = 2)

# kpText(kp, chr=seqlevels(kp$genome), y=0.4, x=0, data.panel = 2, r0=0.2, r1=0, col="#444444", label="30x", cex=0.8, pos=2)
# kpAbline(kp, h=0.4, data.panel = 2, r0=0.2, r1=0, col=data.points.colors[3])






# x<-GRanges(TE_fc_out$annotation$Chr, IRanges(TE_fc_out$annotation$Start, TE_fc_out$annotation$End))

# genomic_fc_out$annotation$Start<-sapply(as.numeric(strsplit(as.character(genomic_fc_out$annotation$Start), ";"), "[[", 1))
# genomic_fc_out$annotation$End<-sapply(as.numeric(strsplit(as.character(genomic_fc_out$annotation$End), ";"), "[[", 1))
# y<-GRanges(genomic_fc_out$annotation$Chr, IRanges(genomic_fc_out$annotation$Start, genomic_fc_out$annotation$End))







