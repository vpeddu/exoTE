library(Rsamtools)
library(ggplot2)
library(rtracklayer)
library(tidyverse)
library(plyr)

args = commandArgs(trailingOnly=TRUE)

#read bedfile
print('reading bed file')
bedfile <-import(args[1], format="bed")

#read bamfile
print('reading bamfile')
bam <- scanBam(args[2])

base = args[3]

df <- data.frame(readname = bedfile$name, bed_info = bedfile$NA.8)
df$biotype<-gsub('\\"','',  df$bed_info)

df$biotype<-gsub('gene_type ', '' , gsub(';','',str_extract(df$bed_info,  'gene_type.(.*?);')))

df$biotype<-factor(df$biotype)

freq_table<-as.data.frame(table(df$readname, df$biotype))
freq_table$Freq<-as.numeric(as.character(freq_table$Freq))
freq_table$Var1<-as.character(freq_table$Var1)

aggregated<- freq_table %>% group_by(Var1) %>%  
  mutate(max_score = max(Freq)) %>% 
  ungroup() %>% 
  filter(Freq==max_score)

aggregated<-aggregated[(!duplicated(aggregated$Var1)),]
aggregated$max_score<-NULL

plot_df <- aggregate(aggregated$Freq, by=list(Category=aggregated$Var2), FUN=sum)
  
save.image(file=paste0(base,"biotyping.Rdata"))

biotype_plot<-ggplot(plot_df, aes(x = Category, y = x, fill = Category)) + 
  geom_bar(stat = 'identity') + 
  theme_classic() + 
  #scale_y_log10(limits = c(1,1e7)) +
  scale_y_log10(limits = c(1, round_any(max(plot_df$x), 1e7, f = ceiling)), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        legend.position = 'none') +
  ylab("Frequency") + 
  xlab('Biotype') +
  scale_fill_viridis_d() 
biotype_plot
ggsave(plot = biotype_plot, file = paste0(base,'biotype_plot.pdf'), height = 5, width = 5)

#Pulling read length by biotype
bam_df<-do.call("DataFrame", bam)
bam_indices<-which(bam_df$qname %in% aggregated$Var1)
bam_biotyped<-as.data.frame(bam_df[bam_indices,])
bam_biotyped<-bam_biotyped[!duplicated(bam_biotyped$qname),]

#merging biotyping and read lengths
biotyped_rl<-full_join(bam_biotyped, aggregated, by = c("qname" = "Var1"))

biotype_rl_plot<-ggplot(biotyped_rl, aes(x = Var2, y = qwidth)) + 
  geom_boxplot() +
  theme_classic() + 
  #scale_y_log10(limits = c(1,1e7)) +
  scale_y_log10(limits = c(10,10000), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        legend.position = 'none') +
  ylab("Read length") + 
  xlab('Biotype') +
  scale_fill_viridis_d() 
biotype_rl_plot

ggsave(plot = biotype_rl_plot, file = paste0(base,'biotype_rl_plot.pdf'), height = 5, width = 5)
