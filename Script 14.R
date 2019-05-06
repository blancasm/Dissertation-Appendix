
######## SCRIPT 14 - Raj et al. ########

# CD4 data #
eu_cd4 <- read.table(file = 'CIS ONLY/Raj/eu_cd4T.tsv', sep = '\t', header = TRUE)
ea_cd4 <- read.table(file = 'CIS ONLY/Raj/ea_cd4T.tsv', sep = '\t', header = TRUE)
aa_cd4 <- read.table(file = 'CIS ONLY/Raj/aa_cd4T.tsv', sep = '\t', header = TRUE)
cd4 <- rbind(eu_cd4,ea_cd4,aa_cd4)
cd4$condition <- "CD4"
cd4$fdr <- 0.05

# monocyte data #
eu_mono <- read.table(file = 'CIS ONLY/Raj/eu_monocytes.tsv', sep = '\t', header = TRUE)
ea_mono <- read.table(file = 'CIS ONLY/Raj/ea_monocytes.tsv', sep = '\t', header = TRUE)
aa_mono <- read.table(file = 'CIS ONLY/Raj/aa_monocytes.tsv', sep = '\t', header = TRUE)
mono <- rbind(eu_mono,ea_mono,aa_mono)
mono$condition <- "mono"
mono$fdr <- 0.05

# Join and format
stack <- rbind(cd4, mono)
stack$type <- "cis"
study_data <- stack[,c(1,6,2,3,9,10,12,11,13)]
colnames(study_data) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","condition","type")  
study_data$study <- "raj_2014"

# Save
write.csv(study_data, "raj_data_stacked.csv", row.names = FALSE)
