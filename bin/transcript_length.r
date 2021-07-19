library(ggplot2)
library(ggpmisc)
library(viridisLite)
library(tidyverse)

files <- list.files(getwd(), pattern = '*flag260.counts.txt', full.names = TRUE)
for(i in 1:length(files)){ 
  if( i == 1){ 
    fdf <- read.csv(files[i], sep = '\t',header = FALSE)
    fdf$sample<-strsplit(basename(files[i]), '[.]')[[1]][1]
  }
  else{ 
    temp<-read.csv(files[i], sep = '\t', header = FALSE)
    temp$sample<-strsplit(basename(files[i]), '[.]')[[1]][1]
    fdf<-rbind(fdf, temp)
  }
}

colnames(fdf)<-c('size','sample')

transcript_lengths<-ggplot(fdf,aes(x = size, color = sample)) + 
  ggtitle('Length of uniquely mapped reads (flag: 260)') +
  scale_color_viridis_d() + 
  geom_freqpoly(binwidth = 10) + 
  scale_y_log10(limits = c(1,1e5)) + 
  #xlim(0,10000) + 
  scale_x_continuous(breaks = seq(0, 10000, by = 1000), limits = c(0,10000)) +
  theme_classic() + 
  geom_vline(xintercept = seq(0, 10000, by = 500),linetype="dotted", color = 'grey') +
  xlab('read length') + 
  ylab('frequency') 

transcript_lengths

ggsave(plot = transcript_lengths, file = 'unnanotated_uniquely_mapped_transcript_lengths.pdf', height = 5, width = 10)

pb <- ggplot_build(transcript_lengths)
labeled_peaks<-transcript_lengths + stat_peaks(
  data = pb[['data']][[1]], # take a look at this object
  aes(x = x, y = count),
  colour = "red",
  size = 3,
  strict = TRUE,
  geom = "text",
  ignore_threshold = 0.3,
  span = 20,
  hjust = -0.1
)
labeled_peaks
ggsave(plot = pb , file = 'peaks_labeled_uniquely_mapped_transcript_lengths.pdf', height = 5, width = 10)




files <- list.files(getwd(), pattern = '*flag2048.counts.txt', full.names = TRUE)
for(i in 1:length(files)){ 
  if( i == 1){ 
    fdf <- read.csv(files[i], sep = '\t',header = FALSE)
    fdf$sample<-strsplit(basename(files[i]), '[.]')[[1]][1]
  }
  else{ 
    temp<-read.csv(files[i], sep = '\t', header = FALSE)
    temp$sample<-strsplit(basename(files[i]), '[.]')[[1]][1]
    fdf<-rbind(fdf, temp)
  }
}

colnames(fdf)<-c('size','sample')

transcript_lengths<-ggplot(fdf,aes(x = size, color = sample)) + 
  ggtitle('Length of supplementary reads (flag: 2048)') +
  scale_color_viridis_d() + 
  geom_freqpoly(binwidth = 10) + 
  scale_y_log10(limits = c(1,1e5)) + 
  #xlim(0,10000) + 
  scale_x_continuous(breaks = seq(0, 10000, by = 1000), limits = c(0,10000)) +
  theme_classic() + 
  geom_vline(xintercept = seq(0, 10000, by = 500),linetype="dotted", color = 'grey') +
  xlab('read length') + 
  ylab('frequency') 

transcript_lengths

ggsave(plot = transcript_lengths, file = 'unnanotated_supplementary_transcript_lengths.pdf', height = 5, width = 10)

pb <- ggplot_build(transcript_lengths)
labeled_peaks<-transcript_lengths + stat_peaks(
  data = pb[['data']][[1]], # take a look at this object
  aes(x = x, y = count),
  colour = "red",
  size = 3,
  strict = TRUE,
  geom = "text",
  ignore_threshold = 0.3,
  span = 20,
  hjust = -0.1
)
labeled_peaks
ggsave(plot = pb , file = 'peaks_labeled_supplementary_transcript_lengths.pdf', height = 5, width = 10)
