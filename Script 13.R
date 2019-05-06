
####### SCRIPT 13 - Nedelec et al. ########

# ALL DATA
cis_data4 <- read.csv('CIS ONLY/Nedelec/cis_eqtl.csv', header = TRUE, skip = 3)
cis_data4$test_beta_Z <- NA
cis_data4$fdr <- 0.05
cis_data4$position <- NA
cis_data4$probe_id <- NA

# Baseline
baseline <- cis_data4[,c(3,27,2,28,25,9,26)]
baseline$condition <- "none"
colnames(baseline) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","condition")

# Listeria
listeria <- cis_data4[,c(4,27,2,28,25,10,26)]
listeria$condition <- "listeria"
colnames(listeria) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","condition")

# Salmonella
salmonella <- cis_data4[,c(5,27,2,28,25,11,26)]
salmonella$condition <- "salmonella"
colnames(salmonella) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","condition")

# Stack and save
cis <- rbind(baseline,listeria,salmonella)
cis$type <- "cis"
cis$study <- "nedelec_2016"
write.csv(cis, "nedelec_data_stacked.csv", row.names = FALSE)
