library("BiocManager") #cite
library("pasilla") #cite?
library(readr)
library("DESeq2") #cite
library("apeglm") #cite
library("RColorBrewer")
library('pheatmap')
library('ggplot2')
library('patchwork')
library('openxlsx')
library('readxl')
#import mature_counts
all <- read.delim("deseq_all.tsv",sep = '\t')
all <-all[-c(8702),]
all <-all[-c(4760),]
rownames(all) <- all$Geneid
all <- all[-1]
## ^^duplicate names exist
smrna <- read.delim("deseq_all_smrna.tsv",sep = '\t',row.names = 1)
mirna <-read.delim("deseq_mirna.tsv",sep = '\t',row.names = 1)
pirna <-read.delim("deseq_pirna.tsv",sep = '\t',row.names = 1)
#create cohorts
comp_25 <- c('CEL1', 'CEL2', 'MED1', 'MED2')
comp_26 <- c('CEL3', 'CEL4', 'MED3','MED4')
comp_27 <- c('CEL5','CEL6','MED5','MED6')
comp_28 <- c('CEL7','CEL8','MED7','MED8')

deseq_inital <- function(dataset){
  df_25 <- dataset[comp_25]
  df_26 <- dataset[comp_26]
  df_27 <- dataset[comp_27]
  df_28 <- dataset[comp_28]
  
  cts_25 <- as.matrix(df_25)
  cts_26 <- as.matrix(df_26)
  cts_27 <- as.matrix(df_27)
  cts_28 <- as.matrix(df_28)
  
  metadata <- read.csv('metadata.csv', row.names = 1)
  
  md_25 <- metadata[rownames(metadata) %in% comp_25, ]
  md_26 <- metadata[rownames(metadata) %in% comp_26, ]
  md_27 <- metadata[rownames(metadata) %in% comp_27, ]
  md_28 <- metadata[rownames(metadata) %in% comp_28, ]
  
  md_25 <- md_25[,c('media', 'cell', 'fbs_exosome', 'type')]
  md_25$media <- factor(md_25$media)
  md_25$cell <- factor(md_25$cell)
  md_25$fbs_exosome <- factor(md_25$fbs_exosome)
  md_25$type <- factor(md_25$type)
  
  md_26 <- md_26[,c('media', 'cell', 'fbs_exosome', 'type')]
  md_26$media <- factor(md_26$media)
  md_26$cell <- factor(md_26$cell)
  md_26$fbs_exosome <- factor(md_26$fbs_exosome)
  md_26$type <- factor(md_26$type)
  
  md_27 <- md_27[,c('media', 'cell', 'fbs_exosome', 'type')]
  md_27$media <- factor(md_27$media)
  md_27$cell <- factor(md_27$cell)
  md_27$fbs_exosome <- factor(md_27$fbs_exosome)
  md_27$type <- factor(md_27$type)
  
  md_28 <- md_28[,c('media', 'cell', 'fbs_exosome', 'type')]
  md_28$media <- factor(md_28$media)
  md_28$cell <- factor(md_28$cell)
  md_28$fbs_exosome <- factor(md_28$fbs_exosome)
  md_28$type <- factor(md_28$type)
  
  num_25 <- matrix(as.numeric(cts_25), ncol = ncol(cts_25))
  rownames(num_25) <- rownames(cts_25)
  colnames(num_25) <- colnames(cts_25)
  
  num_26 <- matrix(as.numeric(cts_26), ncol = ncol(cts_26))
  rownames(num_26) <- rownames(cts_26)
  colnames(num_26) <- colnames(cts_26)
  
  num_27 <- matrix(as.numeric(cts_27), ncol = ncol(cts_27))
  rownames(num_27) <- rownames(cts_27)
  colnames(num_27) <- colnames(cts_27)
  
  num_28 <- matrix(as.numeric(cts_28), ncol = ncol(cts_28))
  rownames(num_28) <- rownames(cts_28)
  colnames(num_28) <- colnames(cts_28)
  
  #
  dds_25 <- DESeqDataSetFromMatrix(countData = num_25,
                                colData = md_25,
                                design = ~ type)
  keep_25 <- rowSums(counts(dds_25)) >= 10
  dds_25 <- dds_25[keep_25,]

  dds_25 <- DESeq(dds_25)
  
  
  #
  dds_26 <- DESeqDataSetFromMatrix(countData = num_26,
                                   colData = md_26,
                                   design = ~ type)
  keep_26 <- rowSums(counts(dds_26)) >= 10
  dds_26 <- dds_26[keep_26,]
  
  dds_26 <- DESeq(dds_26)
  
  #
  dds_27 <- DESeqDataSetFromMatrix(countData = num_27,
                                   colData = md_27,
                                   design = ~ type)
  keep_27 <- rowSums(counts(dds_27)) >= 10
  dds_27 <- dds_27[keep_27,]
  
  dds_27 <- DESeq(dds_27)
  
  #
  dds_28 <- DESeqDataSetFromMatrix(countData = num_28,
                                   colData = md_28,
                                   design = ~ type)
  keep_28 <- rowSums(counts(dds_28)) >= 10
  dds_28 <- dds_28[keep_28,]
  
  dds_28 <- DESeq(dds_28)
  
  
  return_dds <- list(dds_25, dds_26, dds_27, dds_28)
  
  return(return_dds)
}

all_dds <- deseq_inital(all)
mirna_dds <- deseq_inital(mirna)
pirna_dds <- deseq_inital(pirna)
smrna_dds <- deseq_inital(smrna)


##
all_res <- list()
n = 1
for (dds in all_dds){
  all_res[[n]] <- results(dds)
  n = n + 1
}

mirna_res <- list()
n = 1
for (dds in mirna_dds){
  mirna_res[[n]] <- results(dds)
  n = n + 1
}

pirna_res <- list()
n = 1
for (dds in pirna_dds){
  pirna_res[[n]] <- results(dds)
  n = n + 1
}

smrna_res <- list()
n = 1
for (dds in smrna_dds){
  smrna_res[[n]] <- results(dds)
  n = n + 1
}
##

all_rld <- list()
n = 1
for (dds in all_dds){
  all_rld[[n]] <- rlog(dds, blind = FALSE, fitType = 'local')
  n = n + 1
}

mirna_rld <- list()
n = 1
for (dds in mirna_dds){
  mirna_rld[[n]] <- rlog(dds, blind = FALSE, fitType = 'local')
  n = n + 1
}

pirna_rld <- list()
n = 1
for (dds in pirna_dds){
  pirna_rld[[n]] <- rlog(dds, blind = FALSE, fitType = 'local')
  n = n + 1
}

smrna_rld <- list()
n = 1
for (dds in smrna_dds){
  smrna_rld[[n]] <- rlog(dds, blind = FALSE, fitType = 'local')
  n = n + 1
}
rm(dds)

##
#key for numbers
cell_key <- list('SK-N-SH','IMR32','T47D','HCC1143')
##
all_ma <- list()
n = 1
for (res in all_res){
  pdf(file = paste(cell_key[n], 'all.pdf', sep=''))
  plotMA(res, ylim=c(-5.5,5.5), main = paste(cell_key[n], '- Whole Genome Annotation, Media vs. Cell'))
  n = n + 1
  dev.off()
}

mirna_ma <- list()
n = 1
for (res in mirna_res){
  pdf(file = paste(cell_key[n], 'mirna.pdf', sep='_'))
  plotMA(res, ylim=c(-5.5,5.5), main = paste(cell_key[n], '- miRNA specific Annotation, Media vs. Cell'))
  mirna_ma[[n]] <- recordPlot()
  n = n + 1
  dev.off()
}

pirna_ma <- list()
n = 1
for (res in pirna_res){
  pdf(file = paste(cell_key[n], 'pirna.pdf', sep='_'))
  plotMA(res, ylim=c(-5.5,5.5), main = paste(cell_key[n], '- piRNA specific Annotation, Media vs. Cell'))
  pirna_ma[[n]] <- recordPlot()
  n = n + 1
  dev.off()
}

smrna_ma <- list()
n = 1
for (res in smrna_res){
  pdf(file = paste(cell_key[n], 'smrna.pdf', sep='_'))
  plotMA(res, ylim=c(-5.5,5.5), main = paste(cell_key[n], '- smRNA specific Annotation, Media vs. Cell'))
  smrna_ma[[n]] <- recordPlot(n)
  n = n + 1
  dev.off()
}

## PCA

cohort_list <- list(all_rld, mirna_rld, pirna_rld, smrna_rld)
cohort_key <- list('- Whole Genome Annotation',
                   '- miRNA specific Annotation',
                   '- piRNA specific Annotation',
                   '- smRNA specific Annotation')

pca_list <- list()
p = 1
q = 1
n = 1
for (cohort in cohort_list){
  q = 1
  for (rld in cohort){
    pca_list[[p]] <- plotPCA(rld, intgroup=c("fbs_exosome", "type")) + 
      ggtitle(paste('PCA',cell_key[q], cohort_key[n])) + 
      coord_fixed((ratio = 2.5)) +
      ylim(-15,15) + xlim(-50, 50)
    print(q)
    print(n)
    q = q + 1
    p = p + 1
  }
  n = n + 1
}

PCA_SK <- (pca_list[[1]]|pca_list[[2]])/(pca_list[[3]]|pca_list[[4]])
PCA_SK <- PCA_SK + plot_annotation(title = 'Principle Component Analysis for Breast & Neuroblastoma Cell Lines',
                                   subtitle = 'Based upon smRNA-seq alignment with full human genome',
                                   caption = 'Factor 1: plus/minus; presence or absence of FBS exosomes \n
                                   Factor 2: Sample from media or cell')
ggsave('PCA_SK.pdf', plot = PCA_SK, width = 10, height = 6)
PCA_SK <- (pca_list[[1]]|pca_list[[5]])/(pca_list[[9]]|pca_list[[13]])
PCA_SK <- PCA_SK + plot_annotation(title = 'Principle Component Analysis for Breast & Neuroblastoma Cell Lines',
                                   subtitle = 'Based upon smRNA-seq alignment with full human genome',
                                   caption = 'Factor 1: plus/minus; presence or absence of FBS exosomes \n
                                   Factor 2: Sample from media or cell')
                                  

###
cohort_list_res <- list(all_res, mirna_res, pirna_res, smrna_res)
cohort_key_short <- list('- Whole',
                   '- miRNA',
                   '- piRNA',
                   '- smRNA')

ord_list <- list()
p = 1
q = 1
n = 1
for (cohort in cohort_list_res){
  q = 1
  for (res in cohort){
    ord_list[[paste(cell_key[q], cohort_key_short[n])]] <- data.frame(res[order(res$pvalue),])
    q = q + 1
    p = p + 1
  }
  n = n + 1
}

write.xlsx(ord_list, 'DEGs.xlsx', row.names=TRUE)
