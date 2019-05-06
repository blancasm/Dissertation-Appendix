
######## SCRIPT 12 - Manry et al. #######

# Data
cis_stim <- read.csv('CIS ONLY/Manry/cis.csv', header = TRUE, skip = 0)
cis_nonstim <- read.csv('CIS ONLY/Manry/cis_nonstimulated.csv', header = TRUE, skip = 0)
cis_stim$condition <- "M_leprae"
cis_nonstim$condition <- "baseline"

# Stack
cis <- rbind(cis_nonstim, cis_stim)

# Standardise
cis$fdr <- 0.05
cis$probe_id <- NA
cis <- cis[,c(4,6,2,11,7,8,10,9)]
colnames(cis) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr", "condition")  
cis$type <- "cis"
cis$study <- "manry_2017"

# Save
write.csv(cis, "manry_data_stacked.csv", row.names = FALSE)
