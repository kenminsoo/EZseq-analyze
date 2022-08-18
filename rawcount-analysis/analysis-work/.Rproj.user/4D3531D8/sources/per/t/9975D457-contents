library(readr)
library(tidyverse)
library(dplyr)
library('openxlsx')
library('readxl')
library('ggplot2')

sheets <- openxlsx::getSheetNames('tpm_all_sd.xlsx')
data_frame <- lapply(sheets, openxlsx::read.xlsx, xlsxFile='tpm_all_sd.xlsx')

# assigning names to data frame
names(data_frame) <- sheets

adj_list10 = list(3,4)
adj_list12 = list(5,6)

adj3 <- merge(data_frame[[3]], data_frame[[1]], by = 'Geneid')

adj3$`output/aligned/trimmed_merged_Med1-394078718.sam` <- adj3$`output/aligned/trimmed_merged_Med1-394078718.sam` - adj3$`output/aligned/trimmed_merged_Med10-394075788.sam`
adj3[3] <- adj3[3] - adj3[6]

adj3_filt <- adj3 %>%
  filter_all( all_vars(. > 0))

adj3_filt$mean.x <- rowMeans(adj3_filt[2:3])

FBS_minus <- function(adj_num, dataset){  #10 = 6  #12 = 7
  print(data_frame)
  return_df <- merge(data_frame[[dataset]], data_frame[[1]], by = 'Geneid')
  return_df[4] <- return_df[4] - return_df[adj_num]
  return_df[3] <- return_df[3] - return_df[adj_num]
  return_df <- return_df %>%
    filter_all(all_vars(. > 0))
  return_df$mean.x <- rowMeans(return_df[2:3])
  return(return_df[1:5])
}

#this is terrible but i would like to finish
extract_comparison_tables <- function(output_file){
  
  
  #make adjust
  adj3 <- FBS_minus(6, 3)
  adj4 <- FBS_minus(6, 4)
  adj5 <- FBS_minus(7, 5)
  adj6 <- FBS_minus(7, 6)
  
  return_list <- list()
  #create comparison table
  return_list[['comp_25']] <- merge(adj3, data_frame[[7]], by = 'Geneid')
  return_list[['comp_26']] <- merge(adj4, data_frame[[8]], by = 'Geneid')
  return_list[['comp_27']] <- merge(adj5, data_frame[[9]], by = 'Geneid')
  return_list[['comp_28']] <- merge(adj6, data_frame[[10]], by = 'Geneid')

  #filter comparison table
  n = 1
  for (dataset in return_list){
    range_store <- dataset$mean.x - dataset$mean
    print(return_list[[n]])
    print(dataset)
    return_list[[n]] <- cbind(return_list[[n]], range = range_store)
    return_list[[n]] <- return_list[[n]]%>%
      filter(range > 5)
      
    n = n + 1
  }
  
  
  write.xlsx(return_list, output_file)
}

sheets <- openxlsx::getSheetNames('tpm_all_sd.xlsx')
data_frame <- lapply(sheets, openxlsx::read.xlsx, xlsxFile='tpm_all_sd.xlsx')

# assigning names to data frame
names(data_frame) <- sheets

extract_comparison_tables('comp_all.xlsx')
#####
sheets <- openxlsx::getSheetNames('tpm_mirna_sd.xlsx')
data_frame <- lapply(sheets, openxlsx::read.xlsx, xlsxFile='tpm_mirna_sd.xlsx')

# assigning names to data frame
names(data_frame) <- sheets

extract_comparison_tables('comp_mirna.xlsx')
####
sheets <- openxlsx::getSheetNames('tpm_pirna_sd.xlsx')
data_frame <- lapply(sheets, openxlsx::read.xlsx, xlsxFile='tpm_pirna_sd.xlsx')

# assigning names to data frame
names(data_frame) <- sheets

extract_comparison_tables('comp_pirna.xlsx')
####
sheets <- openxlsx::getSheetNames('tpm_smrna_sd.xlsx')
data_frame <- lapply(sheets, openxlsx::read.xlsx, xlsxFile='tpm_smrna_sd.xlsx')

# assigning names to data frame
names(data_frame) <- sheets

extract_comparison_tables('comp_smrna.xlsx')





