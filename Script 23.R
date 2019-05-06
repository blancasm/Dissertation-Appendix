
################ SCRIPT 23 - Figure 5 determination and plot ####################

########################### corr coeff plot data ###############################
truth_results_over20 <- read.csv("Fuzzy/truth_fuzzy_with_mod_over20_MEANS.csv")
truth_results_under20 <- read.csv("Fuzzy/truth_fuzzy_with_mod_under20_MEANS.csv")
truth_results <- rbind(truth_results_over20, truth_results_under20)
truthresults_2col <- truth_results[,c("inter_id", "inter_size","df_means")]
truthresults_2col <- unique(truthresults_2col)


########################### Overlaying for significance ####################
randomco <- read.csv("percentile_random_10000_50_ protect _.csv", header = T)
randomco %>%
  filter(mod_size <= 71) -> randomco_reduced 
randomco_reduced$sd_97 <- (randomco_reduced$X97.72499. - randomco_reduced$X50.)/2
randomco_reduced$sd_84 <- (randomco_reduced$X84.13447. - randomco_reduced$X50.)/1
randomco_reduced$sd_99 <- (randomco_reduced$X99.86501. - randomco_reduced$X50.)/3
randomco_reduced$sd_15 <- (randomco_reduced$X50. - randomco_reduced$X15.86553.)/1
randomco_reduced$sd_2 <- (randomco_reduced$X50. - randomco_reduced$X2.275013.)/2
randomco_reduced$sd <- (randomco_reduced$sd_84 + randomco_reduced$sd_97 + randomco_reduced$sd_99 + randomco_reduced$sd_15 + randomco_reduced$sd_2)/5

######## FIND THE NUMBER OF TRUTH OUTDEGREES BEYOND CERTAIN THRESHOLDS ###########
new_truth_matrix <- matrix(ncol = 8)
colnames(new_truth_matrix) <- c(colnames(truth_results), "significance")
for (i in unique(truth_results$inter_id)){
  truth_results %>%
    filter(inter_id == i) -> hub
  size <- hub$inter_size
  truth <- hub$df_means
  if (size <= 51){
    threshold <- randomco[(size-1),7]
    if (threshold <= truth){
      hub$significance <- TRUE
    } else if (threshold > truth){
      hub$significance <- FALSE
    } else{
      hub$significance <- NA
    }
    new_truth_matrix <- rbind(new_truth_matrix, hub)
  } else if (size == 56){
    if (threshold <= truth){
      hub$significance <- TRUE
    } else if (threshold > truth){
      hub$significance <- FALSE
    } else{
      hub$significance <- NA
    }
    new_truth_matrix <- rbind(new_truth_matrix, hub)
  } else{
    hub$significance <- NA
    new_truth_matrix <- rbind(new_truth_matrix, hub)
  }
}
new_truth_matrix %>%
  filter(significance == TRUE) -> significant_truth_matrix #1,180 individual interactions
significant_cluster_list <- unique(significant_truth_matrix$inter_id) #83 sig
significant_interactome_info <- significant_truth_matrix[, c("inter_id", "target", "source", "inter_size", "df_means")]
write.csv(significant_interactome_info, "sig_inter_indiv_interactions.csv", row.names = F) # IMP FOR P-VAL #

# FIG. 5A
truth_plot <- ggplot(new_truth_matrix, aes(df_means, inter_size, fill = "significance")) +
  geom_point(color = "grey", aes(color = significance)) + xlab("Aveage correlation coefficient") + ylab("Interactome size (number of genes)") + ggtitle("Outdegree Correlation") +
  geom_point(data = randomco_reduced, aes(X97.72499., mod_size), color = "palegreen3") + 
  geom_point(data = randomco_reduced, aes(X84.13447., mod_size), color = "steelblue2") +
  geom_point(data = randomco_reduced, aes(X99.86501., mod_size), color = "sienna2") +
  theme_light(base_size = 10) + theme(legend.position = c(0.75, 0.9)) + 
  scale_fill_discrete(labels = "p-value < FDR threshold")


##### FDR THRESHOLD ####
significant_interactomes_83 <- read.csv("significant_interactome_truth.csv", header = T)
significant_interactomes_83$Z <- (significant_interactomes_83$df_means - significant_interactomes_83$mean_50)/significant_interactomes_83$sd
significant_interactomes_83$pvalue <- pnorm(-abs(significant_interactomes_83$Z))
write.csv(significant_interactomes_83, "significant_interactome_truth_pvalue.csv", row.names = F)

significant_interactomes_truth_pvalue <- subset(significant_interactomes_83, pvalue <= (0.05/83))
sig_hubs_39 <- as.character(significant_interactomes_truth_pvalue$inter_id)
significant_interactomes_83$signif <- NA
for (i in (1:83)){
  if (significant_interactomes_83[i,2] %in% sig_hubs_39){
    significant_interactomes_83[i,11] <- TRUE
  }
}

# FIG. 5B
new_thrsh_plot_2 <- ggplot(significant_interactomes_83, aes(pvalue, inter_size, colour = signif)) +
  geom_point() + xlab("Interactome p-value") + ylab("Interactome size (number of genes)") + ggtitle("Mean Outdegree Truth Significance") + 
  scale_x_continuous(trans='log2', breaks = c(0.000000001,0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.05) ) + 
  geom_vline(xintercept=(0.05/83), colour = "darkred") + 
  theme_light(base_size = 11) +
  labs(fill = "Significance") + 
  scale_fill_discrete(labels = c("p-value > FDR threshold", "p-value < FDR threshold"))

