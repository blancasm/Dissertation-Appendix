
###### SCRIPT 2 - Quach et al. #######

library(tidyr)

# TRANS 
trans <- read.csv("Quach 2016/trans.csv", header = T, skip = 15)

# down
trans_select_down <- trans[,c(1,2,15,3)]
trans_select_down_separated <- separate_rows(trans_select_down, trans.regulated.genes..down.e, sep = " // ", convert = FALSE)
trans_select_down_separated$test_beta_Z <- -1
colnames(trans_select_down_separated) <- c("snp","position","gene_symb","condition", "test_beta_Z")

# up
trans_select_up <- trans[,c(1,2,14,3)]
trans_select_up_separated <- separate_rows(trans_select_up, trans.regulated.genes..up.e, sep = " // ", convert = FALSE)
trans_select_up_separated$test_beta_Z <- 1
colnames(trans_select_up_separated) <- c("snp","position","gene_symb","condition", "test_beta_Z")

# stack trans 
trans_data_stacked <- rbind(trans_select_down_separated, trans_select_up_separated)
# standardise 
trans_data_stacked$type <- "trans"
trans_data_stacked$fdr <- 0.05
trans_data_stacked$p_value <- NA
trans_data_stacked$probe_id <- NA
trans_data_stacked_ordered <- trans_data_stacked[,c(1,2,3,9,5,8,7,4,6)]


# CIS 
cis <- read.csv("Quach 2016/cis.csv", header = T, skip = 13)
cis$probe_id <- NA
cis$type <- "cis"
cis_select <- cis[, c(1,3,5,28,13,11,12,6,29)]
colnames(cis_select) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","condition","type")

# Stack and save
quach_data_stacked <- rbind(cis_select,trans_data_stacked_ordered)
quach_data_stacked$study <- "quach_2016"
write.csv(quach_data_stacked, "quach_data_stacked.csv", row.names = FALSE)
