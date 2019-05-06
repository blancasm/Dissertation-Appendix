
###### SCRIPT 8 - Schramm et al. #######

# CIS
cis <- read.csv('CIS ONLY/Schramm/schramm_cis.csv', header = TRUE, skip = 0)
cis$position <- NA
cis$fdr <- NA
cis <- cis[,c(1,20,6,5,8,11,21)]
colnames(cis) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr")
cis$condition <- "none"
cis$type <- "cis"


# TRANS
trans <- read.csv('CIS ONLY/Schramm/schramm_trans.csv', header = TRUE, skip = 0)
trans$position <- NA
trans$fdr <- NA
trans <- trans[,c(2,19,5,4,7,10,20)]
colnames(trans) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr")
trans$condition <- "none"
trans$type <- "trans"

# Stack and save
data <- rbind(cis, trans)
data$study <- "schramm_2014"
write.csv(data, "schramm_data_stacked.csv", row.names = FALSE)
