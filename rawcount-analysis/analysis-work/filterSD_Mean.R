library(readr)
library(tidyverse)
library('openxlsx')
library('readxl')
library('ggplot2')

#import
tpm_mirna <- read_delim('tpm_mirna.tsv', delim = "\t", escape_double = FALSE, trim_ws = TRUE)
tpm_smrna <- read_delim('tpm_smrna.tsv', delim = "\t", escape_double = FALSE, trim_ws = TRUE)
tpm_pirna <- read_delim('tpm_pirna.tsv', delim = "\t", escape_double = FALSE, trim_ws = TRUE)
tpm_all <- read_delim('tpm_all.tsv', delim = "\t", escape_double = FALSE, trim_ws = TRUE)

#filter
##Please note, do not compare values between datasets

filt_smrna <- tpm_smrna[rowMeans(tpm_smrna[2:24]) > (mean(unlist(tpm_smrna[2:24])) - 5), ] #17.56 - 5 #may have to expand to include more
#filt_smrna2 <- filt_smrna[rowMeans(filt_smrna[2:24]) > mean(unlist(filt_smrna[2:24])), ]

filt_all <- tpm_all[rowMeans(tpm_all[2:24]) > (mean(unlist(tpm_all[2:24])) - 5), ] #17.492 - 5
#filt_all2 <- filt_all[rowMeans(filt_all[2:24]) > mean(unlist(filt_all[2:24])), ]
write.table(filt_all[1],sep = '\t', "set_genes.tsv",quote=F,row.names=F)


filt_pirna <- tpm_pirna[rowMeans(tpm_pirna[2:24]) > (mean(unlist(tpm_pirna[2:24])) - 20), ] #36.10108 - 20
#filt_pirna2 <- filt_pirna[rowMeans(filt_pirna[2:24]) > mean(unlist(filt_pirna[2:24])), ]


filt_mirna <- tpm_mirna[rowMeans(tpm_mirna[2:24]) > (mean(unlist(tpm_mirna[2:24])) - 150), ] #377 - 150
#filt_mirna2 <- filt_mirna[rowMeans(filt_mirna[2:24]) > mean(unlist(filt_mirna[2:24])), ]


## Now let's create cohorts
fbs_rna <- list(12, 14)
con_Med <- list(13, 15, 16)
SK_Med <- list(11, 17)
IMR_Med <- list(18, 19)
T47D_Med <- list(20, 21)
HCC_Med <- list(22, 23)
HMEC_Med <- list(24) #avoid for now
SK_Cel <- list(2,3)
IMR_Cel <- list(4,5)
T47D_Cel <- list(6,7)
HCC_Cel <- list(8,9)
HMEC_Cel <- list(10) #avoid for now
HMEC <- list(10,24)
cohort_names <- list(fbs_rna,con_Med,SK_Med,IMR_Med,T47D_Med,
                     HCC_Med,SK_Cel,IMR_Cel,T47D_Cel,
                     HCC_Cel,HMEC)

#perhaps a bad practice, but R isn't the most beautiful... sorry.
filt_sd <- function(filt_mean){
  rows <- nrow(filt_mean)
  return_df_list <- list()
  n <- 13
  for (cohort in cohort_names){
    temp_df <- data.frame(matrix(ncol = 0, nrow = rows))
    gene_ids <- filt_mean[1]
    temp_df[ , ncol(temp_df) + 1] <- gene_ids
    g <- 1
    for (num in cohort){
      g <- g + 1
      temp_col <- filt_mean[num]
      temp_df[ , ncol(temp_df) + 1] <- temp_col
    }
    temp_sd <- apply(temp_df[2:g], 1, sd)
    temp_df <- cbind(temp_df, sd = temp_sd)
    temp_mean <- apply(temp_df[2:g], 1, mean)
    temp_df <- cbind(temp_df, mean = temp_mean)
    print(temp_df)
    return_df_list[[n]] <- temp_df
    n <- n + 1
  }
  return(return_df_list)
}

tpm_all_list <- filt_sd(filt_all)
tpm_smrna_list <- filt_sd(filt_smrna)
tpm_mirna_list <- filt_sd(filt_mirna)
tpm_pirna_list <- filt_sd(filt_pirna)

test <- tpm_all_list[[12]][tpm_all_list[[12]]['sd'] < median(unlist(tpm_all_list[[12]]['sd'])), ]
test <- test[rowMeans(test[2:3]) > 8, ]

write.xlsx(tpm_all_list[13:23], "tpm_all_sd.xlsx")
write.xlsx(tpm_smrna_list[13:23], "tpm_smrna_sd.xlsx")
write.xlsx(tpm_mirna_list[13:23], "tpm_mirna_sd.xlsx")
write.xlsx(tpm_pirna_list[13:23], "tpm_pirna_sd.xlsx")

#now filtered mean datasets have been created with SD for info. 
##SD will not be used to filter as it will be more selective for lower expressed. Rather, SD is just a way for us to
##better understand our data. 



