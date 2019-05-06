
####### SCRIPT 17 - identification of discrete interaction hubs #######

library(dplyr)

# load files
edges <- read.csv("edges.csv", header = T)
datExpr0 <- read.csv("Fuzzy/tpm_0.001_log2_merge_dedup_symb_ACT_LTB_PPD_for_eqtl.csv", row.names = "X")

# subset so only genes present in the TST expression data are included
base::subset(edges, edges$source %in% colnames(datExpr0) & edges$target %in% colnames(datExpr0)) -> edges_match

# group by common source and put necessary attributes
unique(edges_match$source) -> source_list # number of unique sources - 1,607
hub_data <-setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("source","interaction","target", "weight_cond_scale","weight_study_scale","count","hub_num", "hub_id"))

## For a minimum of 4 genes and a maximum of 20 genes as hub size 
a = 0
for (i in 1:1607){
  edges_match %>%
    filter(source == source_list[i]) -> interactions
  #iteractions[nrow(interactions) + 1,] = list()
  if (4 <= nrow(interactions)) {
    if (nrow(interactions) < 21) {
      a <- a + 1
      interactions %>%
        select(source, interaction, target, weight_cond_scale, weight_study_scale) %>%
        mutate("count" = n()) %>%
        mutate("hub_num" = a) %>%
        mutate("hub_id" = source_list[i]) -> hub_data_1
      hub_data <- rbind(hub_data, hub_data_1)
    }
  }
}

write.csv(hub_data, "blanca_hub_data_20minus.csv", row.names = FALSE) # save file

## For a minimum of 21 genes as hub size
a = 0
for (i in 1:1607){
  edges_match %>%
    filter(source == source_list[i]) -> interactions
  #iteractions[nrow(interactions) + 1,] = list()
  if (21 <= nrow(interactions)) {
     a <- a + 1
     interactions %>%
       select(source, interaction, target, weight_cond_scale, weight_study_scale) %>%
       mutate("count" = n()) %>%
       mutate("hub_num" = a) %>%
       mutate("hub_id" = source_list[i]) -> hub_data_1
     hub_data <- rbind(hub_data, hub_data_1)
  }
}

write.csv(hub_data, "blanca_hub_data_20plus.csv", row.names = FALSE) # save file


