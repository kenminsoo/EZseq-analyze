#load in packages
library("readr")
library("BiocManager") #cite
library("DESeq2") #cite
library("ggplot2")
library('openxlsx')
library('EnhancedVolcano')

#import all files from the count path
import_files <- list.files(path = 'counts', full.names = TRUE)
file_names <- list.files(path = 'counts')
datasets <- list()
n = 1
for (f in import_files){
  datasets[[file_names[n]]] <- read_delim(f, delim = "\t", 
                        escape_double = FALSE, trim_ws = TRUE, 
                        skip = 1)
  n = n + 1
}

raw_col <- ncol(datasets[[1]])

#make the dataset easier to work with
n = 1
for (d in datasets){
  temp <- d
  datasets[[n]] <- data.frame(d)
  rownames(datasets[[n]]) <- temp$Geneid
  datasets[[file_names[n]]] <- datasets[[file_names[n]]][,7:raw_col]
  n = n + 1
}

#clean
rm(d)
rm(thing)
rm(temp)

#import metadata
md <- read.csv('metadata.csv', row.names = 1)

#change to factor
for (name in colnames(md)){
  md[[name]] <- as.factor(md[[name]])
}

#run analyses and export
n = 1
ord_list <- list()
for (d in datasets){
  cts <- as.matrix(datasets[[n]])
  num <- matrix(as.numeric(cts), ncol = ncol(cts))
  rownames(num) <- rownames(cts)
  colnames(num) <- rownames(md)
  
  dds <- DESeqDataSetFromMatrix(countData = num,
                                colData = md,
                                design = ~cell)
  
  dds <- DESeq(dds)
  
  res <- results(dds)
  
  rlog <- rlog(dds, blind = FALSE, fitType = 'local')
  
  pdf(file = paste(file_names[n], '.pdf', sep=''))
  
  plotMA(res, ylim=c(-5.5,5.5), main = paste('RUN', file_names[n]))
  
  dev.off()
  
  pca <- plotPCA(rlog, intgroup=c("type", "cell")) + plot_annotation(title = file_names[n],
                                                              caption = 'Factor 1: Sample from media or cell \n
                                   Factor 2: Cell type')
  
  ggsave((paste(file_names[[n]],'_pca', '.pdf', sep = '')), plot = pca, device = 'pdf')
  
  vol <- EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  title = file_names[n])
  
  ggsave((paste(file_names[[n]],'_vol', '.pdf', sep = '')), plot = vol, device = 'pdf')
  
  
  
  ord_list[[file_names[n]]] <- data.frame(res[order(res$pvalue),])
  n = n + 1
}

write.xlsx(ord_list, 'DEGs.xlsx', rownames=TRUE)
