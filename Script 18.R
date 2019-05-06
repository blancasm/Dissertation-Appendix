
####### SCRIPT 18 - plotting a frequency distribution of hub sizes #######

library(ggplot2)
library(dplyr)

read.csv("blanca_hub_data_20plus.csv") -> edges_20plus
read.csv("blanca_hub_data_20minus.csv") -> edges_20minus
edges <- rbind(edges_20minus, edges_20plus) #before fuzzy clustering

read.csv("all_interactome_data.csv") -> data_fuzzy #after fuzzy clustering

edges %>% # distinct hubs and their counts
  distinct(hub_id, .keep_all = T) -> distinct_hub
distinct_hub <- distinct_hub[,c("count","hub_id")]

data_fuzzy %>% 
  distinct(inter_id, .keep_all = T) -> distinct_hub_fuzzy
distinct_hub_fuzzy <- distinct_hub_fuzzy[,c("inter_size","inter_id")]

# Plot for before fuzzy clusering
freq_distribution <- ggplot(distinct_hub, aes(count)) +
  geom_histogram(fill = "magenta", binwidth = 10) + xlab("Hub size (number of genes)") +
  ylab("Frequency") + ggtitle("Hub size distribution before fuzzy clustering") +
  theme_light(base_size = 12)

# Plot for after fuzzy clustering
freq_distribution_fuzzy <- ggplot(distinct_hub_fuzzy, aes(inter_size)) +
  geom_histogram(fill = "magenta", binwidth = 5) + xlab("Hub size (number of genes)") +
  ylab("Frequency") + ggtitle("Hub size distribution after fuzzy clustering") +
  theme_light(base_size = 12)

# write csv files to use on Rawgraph website platform
write.csv(distinct_hub, "distinct_hub.csv", row.names = F)
write.csv(distinct_hub_fuzzy, "distinct_hub_fuzzy.csv", row.names = F)


# I saved the plot as a PNG using the RStudio iterface
