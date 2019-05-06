
######## SCRIPT 4 - Fairfax et al. ##########

# CIS
cis <- read.csv("Fairfax 2014/cis_data.csv", header = T)
cis$type = "cis"

# Deconstruct the data by condition
cis_LPS2 <- cis[,c(1,5,2,3,7,11,15,26)]
cis_LPS2$condition <- "LPS2"
colnames(cis_LPS2) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","type","condition")

cis_LPS24 <- cis[c(1,5,2,3,8,12,16,26)]
cis_LPS24$condition <- "LPS24"
colnames(cis_LPS24) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","type","condition")

cis_INF <- cis[c(1,5,2,3,9,13,17,26)]
cis_INF$condition <- "INF"
colnames(cis_INF) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","type","condition")

cis_naive <- cis[c(1,5,2,3,10,14,18,26)]
cis_naive$condition <- "naive"
colnames(cis_naive) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","type","condition")

# Stack the different conditions together
cis_stacked <- rbind(cis_LPS2, cis_LPS24, cis_INF, cis_naive)


# TRANS
trans <- read.csv("Fairfax 2014/trans_data.csv", header = T)
trans$type = "trans"

trans_LPS2 <- trans[,c(1,5,2,3,8,12,16,27)]
trans_LPS2$condition <- "LPS2"
colnames(trans_LPS2) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","type","condition")

trans_LPS24 <- trans[c(1,5,2,3,9,13,17,27)]
trans_LPS24$condition <- "LPS24"
colnames(trans_LPS24) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","type","condition")

trans_INF <- trans[c(1,5,2,3,10,14,18,27)]
trans_INF$condition <- "INF"
colnames(trans_INF) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","type","condition")

trans_naive <- trans[c(1,5,2,3,11,15,19,27)]
trans_naive$condition <- "naive"
colnames(trans_naive) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","type","condition")

trans_stacked <- rbind(trans_LPS2, trans_LPS24, trans_INF, trans_naive)


# Stack cis and trans
fairfax_data <- rbind(cis_stacked, trans_stacked)
fairfax_data$study <- "fairfax_2014"


# Filter by FDR < 0.05
fairfax_data_select <- subset(fairfax_data, fairfax_data$fdr < 0.05 )


# Save standardised table data
write.csv(fairfax_data_select, "fairfax_data_stacked.csv", row.names = FALSE)



