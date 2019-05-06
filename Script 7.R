
###### SCRIPT 7 - Naranbhai et al. #######

# CIS
cis <- read.csv("Naranbhaj 2015/2015cis_eqtl.csv", header = T)

cis$type <- "cis"
cis$position <- NA
cis$condition <- "CD16 Neutrophil"

cis_select <- cis[,c(1,10,8,2,6,4,5,11,9)]

# TRANS 
trans <- read.csv("Naranbhaj 2015/2015trans_eqtl.csv", header = T)

trans$type <- "trans"
trans$position <- NA
trans$condition <- "CD16 Neutrophil"
trans$gene <- NA

trans_select <- trans[,c(1,9,7,11,6,4,5,10,8)]

# Stack and save
naranbhai_data <- rbind(cis_select, trans_select)
colnames(naranbhai_data) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","condition","type")
naranbhai_data$study <- "naranbhai_2015"

write.csv(naranbhai_data, "naranbhai_data_stacked.csv", row.names = FALSE)

