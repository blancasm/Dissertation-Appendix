
####### SCRIPT 11 - Kim-Hellmuth et al. ########

cis <- read.csv('CIS ONLY/Kim-Hellmuth/cis.csv', header = TRUE, skip = 0)
cis$test_beta_Z <- NA

# Control
cis_control <- cis[,c(1,5,3,2,28,9)]
colnames(cis_control) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value")
cis_control$fdr <- NA
cis_control$condition <- "none"

# LPS 90 min
cis_lps90 <- cis[,c(1,5,3,2,28,10,22)]
colnames(cis_lps90) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr")
cis_lps90$condition <- "LPS_90min"

# LPS 6 hr
cis_lps6 <- cis[,c(1,5,3,2,28,11,23)]
colnames(cis_lps6) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr")
cis_lps6$condition <- "LPS_6hr"

# RNA 90 min
cis_rna90 <- cis[,c(1,5,3,2,28,12,24)]
colnames(cis_rna90) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr")
cis_rna90$condition <- "RNA_90min"

# RNA 6 hr
cis_rna6 <- cis[,c(1,5,3,2,28,13,25)]
colnames(cis_rna6) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr")
cis_rna6$condition <- "RNA_6hr"

# MDP 90 min
cis_mdp90 <- cis[,c(1,5,3,2,28,14,26)]
colnames(cis_mdp90) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr")
cis_mdp90$condition <- "MDP_90min"

# MDP 6 hr
cis_mdp6 <- cis[,c(1,5,3,2,28,15,27)]
colnames(cis_mdp6) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr")
cis_mdp6$condition <- "MDP_6hr"

# Join and save
cis_stacked <- rbind(cis_control, cis_lps6, cis_lps90, cis_mdp6, cis_mdp90, cis_rna6, cis_rna90)
cis_stacked$type <- "cis"
cis_stacked$study <- "kim-hellmuth_2017" # 17,626 observations
write.csv(cis_stacked, "kimhellmuth_data_stacked.csv", row.names = FALSE) 
