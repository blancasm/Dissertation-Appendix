
####### SCRIPT 5 - Kasela et al. ########

library(tidyr)

# CIS
# CD4 data 
cis1 <- read.csv("Kasela 2017/cis_CD4.csv", header = TRUE, skip = 0)
# CD8 data 
cis2 <- read.csv("Kasela 2017/cis_CD8.csv", header = TRUE, skip = 0)
colnames(cis2)[1] <- "PValue"

# Divide columns with ; into separate columns
cis1 <- separate_rows(cis1, HGNCName, sep = ";", convert = FALSE)
cis2 <- separate_rows(cis2, HGNCName, sep = ";", convert = FALSE)

# Join CD4 and CD8 data 
cis <- rbind(cis1, cis2)
cis$type <- "cis"
cis$condition <- NA
cis$beta <- NA

cis_select <- cis[,c(2,4,13,5,19,1,15,20,18)]
colnames(cis_select) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","condition","type")


# TRANS 
trans1 <- read.csv("Kasela 2017/trans_CD4.csv", skip = 0)
colnames(trans1)[1] <- "PValue"
trans2 <- read.csv("Kasela 2017/trans_CD8.csv", skip = 0)
trans1 <- separate_rows(trans1, HGNCName, sep = ";", convert = FALSE)
trans2 <- separate_rows(trans2, HGNCName, sep = ";", convert = FALSE)

# Join CD4 and CD8 data #
trans <- rbind(trans1, trans2)
trans$type <- "trans"
trans$beta <- NA

# Match columns to those found in cis data #
trans_select <- trans[,c(2,4,13,5,1,22,15,18,21)]
colnames(trans_select) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","condition","type")  

 
## BINDING TRANS AND CIS DATA ##
kasela_data <- rbind(cis_select,trans_select)
kasela_data$study <- "kasela_2015"

write.csv(kasela_data, "kasela_data_stacked.csv", row.names = FALSE)
