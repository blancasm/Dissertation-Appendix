
####### SCRIPT 3 - Davenport et al. paper #######

## Load spreadsheet data

#CIS
cis <- read.csv("Davenport 2016/cis_eqtls.csv", header = T)
cis$X <- NULL
cis$type <- "cis"
cis$condition <- "sepsis"
cis_data <- cis[c(1,3,5,4,6,8,9,11,12)]
colnames(cis_data) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","type","condition")


#TRANS
trans <- read.csv("Davenport 2016/cis_trans_eqtls.csv", header = T)
trans$type <- "trans"
trans$condition <- "sepsis"
trans_data <- trans[c(1,3,7,4,8,10,11,16,17)]
colnames(trans_data) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","type","condition")

#Join and save data
davenport_data <- rbind(cis_data,trans_data)
davenport_data$study <- "davenport_2016"
write.csv(davenport_data, "davenport_data.csv", row.names = FALSE)
