
###### SCRIPT 6 - Kim et al. ########


# LPS stimulation
LPS <- read.csv("Kim 2014/LPS data.csv", header = T)
LPS$condition <- "LPS"
LPS_select <- LPS[c(4,6,1,2,8,7,9,11,10)]
colnames(LPS_select) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","condition","type")

# Baseline
baseline <- read.csv("Kim 2014/baseline data.csv", header = T)
baseline$consition <- "baseline"
baseline_select <- baseline[c(4,6,1,2,8,7,9,11,10)]
colnames(baseline_select) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","condition","type")

# Stack and save
kim_data <- rbind(LPS_select, baseline_select)
kim_data$study <- "kim_2014"
write.csv(kim_data, "kim_data_stacked.csv", row.names = FALSE)
