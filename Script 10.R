
######## SCRIPT 10 - Barreiro et al. paper #########

cisdata <- read.csv('CIS ONLY/Barreiro/barreiro_cis_data.csv', header = TRUE, skip = 2)
cisdata$position <- NA
cisdata$probe_id <- NA
cisdata$test_beta_Z <- NA
cisdata$fdr <- NA

## non-stimulated
cis_unstim <- cisdata[,c(3,18,1,19,20,4,21)]
colnames(cis_unstim) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr")
cis_unstim$condition <- "none"

## stimulated
cis_stim <- cisdata[,c(5,18,1,19,20,6,21)]
colnames(cis_stim) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr")
cis_stim$condition <- "TB"

# join
cis <- rbind(cis_stim, cis_unstim)
cis$type <- "cis"
cis$study <- "barreiro_2012"

# select p-value of less than 0.05
library(dplyr)
cis %>%
  filter(p_value <= 0.05) -> cis_filtered # filters from 23,926 to 19,213

write.csv(cis_filtered, "barreiro_data_stacked.csv", row.names = FALSE)
