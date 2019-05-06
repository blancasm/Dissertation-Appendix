
######## SCRIPT 15 - Stacking eQTL data from all cis-only and cis and trans studies ########

# Load cis and trans study data

naranbhai_data <- read.csv("naranbhai_data_stacked.csv", header = T)
kasela_data <- read.csv("kasela_data_stacked.csv", header = T)
fairfax_data <- read.csv("fairfax_data_stacked.csv", header = T)
kim_data <- read.csv("kim_data_stacked.csv", header = T)
quach_data <- read.csv("quach_data_stacked.csv", header = T)
lee_data <- read.csv("lee_data.csv", header = T)
davenport_data <- read.csv("davenport_data.csv", header = T)
schramm_data <- read.csv("schramm_data_stacked.csv", header = T)

# Load cis-only study data

barreiro_data <- read.csv("barreiro_data_stacked.csv", header = T)
raj_data <- read.csv("raj_data_stacked.csv", header = T)
manry_data <- read.csv("manry_data_stacked.csv", header = T)
kimhellmuth_data <- read.csv("kimhellmuth_data_stacked.csv", header = T)
nedelec_data <- read.csv("nedelec_data_stacked.csv", header = T)

# Stack all study data
stacked_data <- rbind(naranbhai_data, kasela_data, fairfax_data, kim_data, quach_data, lee_data, davenport_data, schramm_data, barreiro_data, nedelec_data, manry_data, kimhellmuth_data, raj_data)
# 2,017,139 individual interactions
  
cis_data_stacked <- stacked_data[stacked_data$type == "cis",] #2,008,156 cis interactions
trans_data_stacked <- stacked_data[stacked_data$type == "trans",] #32,143 trans interactions

cis_data_stacked$type = NULL
trans_data_stacked$type = NULL

# Save
write.csv(cis_data_stacked, "MAPPING PREPARATION/all_studies_cis.csv", row.names = T)
write.csv(trans_data_stacked, "MAPPING PREPARATION/all_studies_trans.csv", row.names = T)
