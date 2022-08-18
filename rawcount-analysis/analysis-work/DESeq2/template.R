library("BiocManager") #cite
library("pasilla") #cite?
library(readr)
library("DESeq2") #cite
library("apeglm") #cite
library("RColorBrewer")
library('pheatmap')
library('ggplot2')
#import mature_counts
mature_counts <- read.csv("mature_counts.csv", row.names = 1)
#SUBSET FOR FIRST COMPARISON

myvars <- c('merged_Med10','merged_Med12','merged_Med6.','merged_Med8.','merged_Med4.','merged_Med2.')
mact_nb <- mature_counts[myvars]

cts <- as.matrix(mact_nb)
#import coldata

mature_coldata <- read.csv("mature_coldata5.csv", row.names = 1)

rownames(mature_coldata) <- sub("-", ".", rownames(mature_coldata))

mact_nb_coldata <- mature_coldata[rownames(mature_coldata) %in% myvars, ]

mact_nb_coldata <- mact_nb_coldata[,c("condition","type")]
mact_nb_coldata$condition <- factor(mact_nb_coldata$condition)
mact_nb_coldata$type <- factor(mact_nb_coldata$type)

#check order
head(cts,2)
mact_nb_coldata
#check names
all(rownames(mact_nb_coldata) == colnames(mact_nb))
#if true, proceed

#convert numeric
mact_nb_num <- matrix(as.numeric(cts), ncol = ncol(cts))
rownames(mact_nb_num) <- rownames(cts)
colnames(mact_nb_num) <- colnames(cts)



dds <- DESeqDataSetFromMatrix(countData = mact_nb_num,
                              colData = mact_nb_coldata,
                              design = ~ condition)

#prefiltering - removed ~1600 mirnas
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

#This will be for FBS-NB vs NB
#deseq analysis
dds <- DESeq(dds)
res <- results(dds)
res
resultsNames(dds)

resOrdered <- res[order(res$pvalue),]
summary(res)
plotMA(res, ylim=c(-5.5,5.5), main = "MA Plot for FBS.pos_vs_FBS.neg")

#significant genes found--conduct further probing
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
colnames(sampleDistMatrix) <- NULL
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, main = "Sample to Sample Distances for FBS.pos_vs_FBS.neg")

plotPCA(rld, intgroup=c("condition", "type")) + ggtitle("Principal Component Analysis of FBS.pos_vs_FBS.neg")

write.csv(resOrdered, "FBS.pos_vs_FBS.neg.csv")
