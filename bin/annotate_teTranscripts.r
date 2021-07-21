library(tidyverse)
library(biomaRt)
library(ensembldb)


#args = commandArgs(trailingOnly=TRUE)

files <- list.files(getwd(), pattern = '*.TEtranscripts.txt', full.names = TRUE)
for(i in 1:length(files)){ 
  if( i == 1){ 
    fdf <- read.csv(files[i], sep = '\t')
    colnames(fdf)[ncol(fdf)]<-strsplit(basename(files[i]), '[.]')[[1]][1]
  }
  else{ 
    temp<-read.csv(files[i], sep = '\t')
    fdf<-cbind(fdf, temp[,2])
    colnames(fdf)[ncol(fdf)]<-strsplit(basename(files[i]), '[.]')[[1]][1]
    }
  }

mapped_norm<-fdf
for(i in 2:ncol(fdf)){ 
  mapped_norm[,i]<-1e6 * fdf[,i] / fdf[(nrow(fdf) - 1),i]
}

mart <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart)
gene_lookup<-getBM(attributes = c('external_gene_name', "ensembl_gene_id_version",'start_position','end_position','chromosome_name'),
                   mart = ensembl)
gene_lookup$chromosome_name<-paste0('chr',gene_lookup$chromosome_name)
gene_lookup$chromosome_name[gene_lookup$chromosome_name == 'chrMT'] <- 'chrM'
#colnames(gene_lookup)[2]<-'gene.TE'

genomic_merged <- dplyr::full_join(mapped_norm,gene_lookup, by = c('gene.TE' = 'ensembl_gene_id_version'))

write.csv(genomic_merged, file = 'annotated_CPM.csv', row.names = FALSE)
