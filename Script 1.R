
########## SCRIPT 1 - TST co-correlation matrix #########

# Load data
tst_metadata <- read.csv("meta_data_raw_salines_06.csv", header = T)
tst_matrix <- read.csv("tpm_0.001_log2_merge_dedup_symb.csv", header = T, row.names = "X")


# Order the matrix and metadata
tst_metadata <- tst_metadata[order(tst_metadata$X),]
tst_matrix <- tst_matrix[, order(colnames(tst_matrix))]

# Computing a correlation matrix of the TST data
write.csv(tst_matrix, "whole_matrix.csv") # saving tst_matrix as a csv file
r_tst <- cor(t(tst_matrix), method = "pearson") 
# r_tst file is too large to save, so in future scripts, load tst_matrix file and then use this code









