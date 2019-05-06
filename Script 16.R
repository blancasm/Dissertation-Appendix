
########## SCRIPT 16 - Generating nodes and edges from the eQTL data ########## 

library(tidyr)
library(dplyr)
library(bpa)

# Import eQTL data
cis_data <- read.csv("all_studies_cis.csv", header = T) # cis eQTL's
trans_data <- read.csv("all_studies_trans.csv", header = T) # trans eQTL's

# Standardize gene names

# make all gene name letters uppercase 
toupper(cis_data$gene_symb) -> cis_data$gene_symb
toupper(trans_data$gene_symb) -> trans_data$gene_symb

# remove NA's or blanks in the snp or gene name
# TRANS
trans_data %>%    # remove NA's
  tidyr::drop_na(snp, gene_symb) -> trans_data
# Replace empty spaces with NA's
for (i in 1:nrow(trans_data)){  
  trans_data$gene_symb[i] %>%
    get_pattern -> pattern
  if (startsWith(pattern, "A") == FALSE) {
      trans_data$gene_symb[i] <- NA
  }
}
trans_data %>%    # remove NA's
  tidyr::drop_na(snp, gene_symb) -> trans_merge

# CIS 
cis_data %>%    # remove NA's
  tidyr::drop_na(snp, gene_symb) -> cis_data
# Replace empty spaces with NA's
for (i in 1:nrow(cis_data)){  
  cis_data$gene_symb[i] %>%
    get_pattern -> pattern
  if (startsWith(pattern, "A") == FALSE) {
    cis_data$gene_symb[i] <- NA
  }
}
cis_data %>%    # remove NA's
  tidyr::drop_na(snp, gene_symb) -> cis_merge 

trans <- trans_merge
cis <- cis_merge

# Identify the SNP's with both cis and trans associations #
SNP_in_both <- intersect(trans$snp, cis$snp) # 4238 

# Filter the data sets so only SNP's in both are present
# CIS
cis_min_snp <- base::subset(cis, cis$snp %in% SNP_in_both) #45605
trans_min_snp <- base::subset(trans, trans$snp %in% SNP_in_both) #27071

cis_min_snp %>%
  group_by(snp, gene_symb) %>%
  summarise(n_conditions = n_distinct(condition[!is.na(condition)])) -> cis_min_snp_cond

cis_min_snp %>%
  group_by(snp, gene_symb) %>%
  summarise(n_study = n_distinct(study[!is.na(study)])) -> cis_min_snp_study
cis_min_snp_cond$n_study <- cis_min_snp_study$n_study

# TRANS
trans_min_snp %>%
  group_by(snp, gene_symb) %>%
  summarise(n_conditions = n_distinct(condition[!is.na(condition)])) -> trans_min_snp_cond

trans_min_snp %>%
  group_by(snp, gene_symb) %>%
  summarise(n_study = n_distinct(study[!is.na(study)])) -> trans_min_snp_study
trans_min_snp_cond$n_study <- trans_min_snp_study$n_study

# Full join
full_join(cis_min_snp_cond, trans_min_snp_cond, by= "snp") -> full_join

# Add weight for study and condition
length(levels(as.factor(cis$condition))) -> n_con
length(levels(as.factor(cis$study))) -> n_study

full_join$weight_cond <- (full_join$n_conditions.x/n_con +full_join$n_conditions.y/n_con) # add weight
full_join$weight_study <- (full_join$n_study.x/n_study +full_join$n_study.y/n_study)

library(scales)
rescale(full_join$weight_cond, to = c(0, 1)) -> full_join$weight_cond_scale   
rescale(full_join$weight_study, to = c(0,1)) -> full_join$weight_study_scale   

# Remove the snp
full_join %>%
  arrange(desc(weight_cond, weight_study)) %>%
  group_by(gene_symb.x, gene_symb.y) %>%
  slice(1:1) -> full_join_unique

# Format edge file for cytoscape 
colnames(full_join_unique) <- 
  gsub(x = colnames(full_join_unique), pattern = "gene_symb.x", replacement = "source")  # cis
colnames(full_join_unique) <- 
  gsub(x = colnames(full_join_unique), pattern = "gene_symb.y", replacement = "target")  # trans
  
# Save
full_join_unique$interaction <- "directed" # required for cytoscape
full_join_unique -> edges # 9.823 unique edges
write.csv(edges, "edges.csv", row.names = FALSE) 



# Format nodes file for cytoscape
as.character(full_join_unique$target) -> full_join_unique$target
as.character(full_join_unique$source) -> full_join_unique$source
c(full_join_unique$source, full_join_unique$target) -> list_nodes
Nodes = unique(list_nodes)
as.data.frame(Nodes) -> nodes
nodes$id <- nodes$Nodes # duplicate ("id": cytoscape, "Id": gephi)
as.character(nodes$Nodes) -> nodes$Nodes

# Add the class
# i.e. cis, trans or both
# 1 = cis only gene
# 2 = trans only gene
# 3 = trans and cis gene
as.vector(full_join_unique$source) -> cis_list
as.vector(full_join_unique$target) -> trans_list

# loop over trans
for (i in 1:nrow(nodes)){
  x <- nodes$Nodes[i]
  nodes$trans[i] <- if_else(x %in% trans_list, 1, 0)
}
# loop over cis
for (i in 1:nrow(nodes)){
  x <- nodes$Nodes[i]
  nodes$cis[i] <- if_else(x %in% cis_list, 1, 0)
}
# add class 
nodes$class <- if_else(nodes$cis ==1 & nodes$trans == 0, 1, 
                       if_else(nodes$trans ==1 & nodes$cis== 0, 2, 3))
                    
write.csv(nodes, "nodes_with_class.csv", row.names = FALSE) 
# 3125 unique nodes
