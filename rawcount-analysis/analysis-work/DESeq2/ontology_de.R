library('openxlsx')
library('readxl')
library('ggplot2')
library('gprofiler2')
library(readr)
library('biomaRt')
library('plotly')
library('dplyr')
#import 
SK <- data_frame[[1]][1]
IMR <- data_frame[[2]][1]
T47 <- data_frame[[3]][1]
HCC <- data_frame[[4]][1]
#select only signifcant
SK_select <- data_frame[[1]][1] %>% filter(data_frame[[1]][7] > .05)
IMR_select <- data_frame[[2]][1] %>% filter(data_frame[[2]][7] > .05)
T47_select <- data_frame[[3]][1] %>% filter(data_frame[[3]][7] > .05)
HCC_select <- data_frame[[4]][1] %>% filter(data_frame[[4]][7] > .05)

sheets <- openxlsx::getSheetNames('Figures_and_Results/DEGs.xlsx')
data_frame <- lapply(sheets, openxlsx::read.xlsx, xlsxFile='Figures_and_Results/DEGs.xlsx')
data_frame <- data_frame[2:5]
# assigning names to data frame
names(data_frame) <- sheets
tpm_all <- read_delim('hsa_all.tsv', delim = '\t', skip = 1)
df_row <- data.frame()
df_row <- rbind(c(tpm_all[1]))

genes <- read_delim('set_genes.tsv', delim = '\t')

##lets analyze what pathways are upregulated in media smRNAs
gost_SK <- gost(query = c(SK_select[1]), organism = 'hsapiens',
                   significant = TRUE,
                   domain_scope = 'custom_annotated',
                   custom_bg = dplyr::pull(SK, X1))

write.xlsx(gost_SK[1], 'ontology_SK.xlsx')

plot1 <- gostplot(gost_SK, interactive = TRUE)
plot1 <- plot1 %>% layout(title = 'Ontology Analysis of SK')

gost_IMR <- gost(query = c(IMR_select[1]), organism = 'hsapiens',
                   significant = TRUE,
                   domain_scope = 'custom_annotated',
                   custom_bg = dplyr::pull(IMR, X1))

write.xlsx(gost_IMR[1], 'ontology_IMR.xlsx')

gost_T47 <- gost(query = c(T47_select[1]), organism = 'hsapiens',
                   significant = TRUE,
                   domain_scope = 'custom_annotated',
                   custom_bg = dplyr::pull(T47, X1))

write.xlsx(gost_T47[1], 'ontology_T47.xlsx')


gost_HCC <- gost(query = c(HCC_select), organism = 'hsapiens',
                   significant = TRUE,
                   domain_scope = 'custom_annotated',
                   custom_bg = dplyr::pull(HCC, X1))

write.xlsx(gost_HCC[1], 'Ontology_HCC.xlsx')









#label it on the sheet
ensembl <- useEnsembl(biomart="genes",dataset='hsapiens_gene_ensembl')

select_genes <- data.frame()
select_genes <- rbind(HCC_select, IMR_select, SK_select, T47_select)
select_genes <- distinct(select_genes)

go_terms <- getBM(mart=ensembl, attributes=c("hgnc_symbol", "uniprot_gn_id", "uniprot_gn_symbol", "go_id", "namespace_1003", "name_1006"),
                  filters="hgnc_symbol", values=c(dplyr::pull(select_genes, X1)))

go_terms$Geneid <- go_terms$hgnc_symbol
go_terms <- go_terms[-1]


HCC_Labeled <- merge(HCC_select, go_terms, by='Geneid')
write.xlsx(HCC_Labeled, 'HCC_Labeled.xlsx')

IMR_Labeled <- merge(IMR_select, go_terms, by='Geneid')
write.xlsx(IMR_Labeled, 'IMR_Labeled.xlsx')

SK_Labeled <- merge(SK_select, go_terms, by='Geneid')
write.xlsx(SK_Labeled, 'SK_Labeled .xlsx')

T47_Labeled <- merge(T47_select, go_terms, by='Geneid')
write.xlsx(T47_Labeled, 'T47_Labeled.xlsx')


