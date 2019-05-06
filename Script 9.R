
####### SCRIPT 9 - Lee et al. ########

# CIS

# Baseline
cis_baseline <- read.csv("Lee 2014/cis-eQTLs.csv", header = T)
cis_baseline$condition = "baseline"

# LPS
cis_LPS <- read.csv("Lee 2014/cis_LPS.csv", header = T)
cis_LPS$condition = "LPS"

# Flu
cis_Flu <- read.csv("Lee 2014/cis_FLU.csv", header = T)
cis_Flu$condition = "FLU"

# IFNb
cis_IFN <- read.csv("Lee 2014/cis_IFNb.csv", header = T)
cis_IFN$condition = "IFNb"

# Join
cis_data <- rbind(cis_baseline, cis_LPS, cis_Flu, cis_IFN)
cis_data$type = "cis"
cis_data$position <- NA
cis_data$probe_id <- NA


# TRANS DATA

# Baseline
trans_baseline <- read.csv("Lee 2014/trans_baseline.csv", header = T)
trans_baseline$condition = "baseline"

# LPS
trans_LPS <- read.csv("Lee 2014/trans_LPS.csv", header = T)
trans_LPS$condition = "LPS"

# Flu
trans_Flu <- read.csv("Lee 2014/trans_FLU.csv", header = T)
trans_Flu$condition = "FLU"

# IFNb
trans_IFN <- read.csv("Lee 2014/trans_IFNb.csv", header = T)
trans_IFN$condition = "IFNb"

# Join
trans_data <- rbind(trans_baseline, trans_LPS, trans_Flu, trans_IFN)
trans_data$type <- "trans"
trans_data$position <- NA
trans_data$probe_id <- NA

# Rename cis
cis_data_renamed <- cis_data[c(2,13,3,14,5,7,8,11,12)]
colnames(cis_data_renamed) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","type","condition")

# Rename trans
trans_data_renamed <- trans_data[c(2,14,3,15,5,7,8,12,13)]
colnames(trans_data_renamed) <- c("snp","position","gene_symb","probe_id","test_beta_Z","p_value","fdr","type","condition")

## Stack and save
lee_data <- rbind(cis_data_renamed, trans_data_renamed)
lee_data$study <- "lee_2014"
write.csv(lee_data, "lee_data.csv", row.names = FALSE)
