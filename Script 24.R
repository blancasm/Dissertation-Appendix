
######## SCRIPT 24 - Figure 6 ##########

library(ggplot2)

# Plot results
random_inter_results <- read.csv("random_single_inter_results.csv", header = T) # load the file where manually loaded query results

freq_distribution_known <- ggplot(random_inter_results, aes(n_reported)) +
  geom_histogram(fill = "purple") + xlab("Number of interactions per 549 gene pairs") +
  ylab("Frequency") + ggtitle("Frequency Distribution of Reported Interactions") +
  theme_light(base_size = 11)
  

freq_distribution_sig <- ggplot(random_inter_results, aes(n_reported)) +
  geom_histogram(fill = "purple") + xlab("Number of interactions per 549 gene pairs") +
  ylab("Frequency") + ggtitle("Frequency Distribution of Reported Significant Interactions") +
  theme_light(base_size = 11)