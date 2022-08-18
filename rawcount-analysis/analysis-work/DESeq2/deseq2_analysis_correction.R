library("BiocManager") #cite
library("pasilla") #cite?
library(readr)
library(tidyverse)
library("DESeq2") #cite
library("apeglm") #cite
library("RColorBrewer")
library('pheatmap')
library('ggplot2')
library('patchwork')
library('openxlsx')
library('readxl')

all <- read.delim("deseq_all.tsv",sep = '\t')
all <-all[-c(8702),]
all <-all[-c(4760),]
rownames(all) <- all$Geneid
all <- all[-1]
## ^^duplicate names exist
smrna <- read.delim("deseq_all_smrna.tsv",sep = '\t',row.names = 1)
mirna <-read.delim("deseq_mirna.tsv",sep = '\t',row.names = 1)
pirna <-read.delim("deseq_pirna.tsv",sep = '\t',row.names = 1)

datasets <- list(all, mirna, pirna, smrna)
corrected <- list()
n = 1 
#correct
for (dataset in datasets){
  corrected[[n]] <- dataset
  corrected[[n]][10] <- dataset[10] - dataset[11]
  corrected[[n]][16] <- dataset[16] - dataset[11]
  corrected[[n]][17] <- dataset[17] - dataset[11]
  corrected[[n]][18] <- dataset[18] - dataset[11]
  
  corrected[[n]][19] <- dataset[10] - dataset[13]
  corrected[[n]][20] <- dataset[20] - dataset[13]
  corrected[[n]][21] <- dataset[21] - dataset[13]
  corrected[[n]][22] <- dataset[21] - dataset[13]
  
  n = n + 1
}

all <- corrected[[1]]
mirna <- corrected[[2]]
pirna <- corrected[[3]]
smrna <- corrected[[4]]


