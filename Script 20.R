
###### SCRIPT 20 - calculating the truth interactome correlation coefficients (individual and mean) ########

library(dplyr)
outdir <- "Fuzzy/"

# Read in the data
datExpr0 = read.csv("Fuzzy/tpm_0.001_log2_merge_dedup_symb_ACT_LTB_PPD_for_eqtl.csv", row.names = "X")
colnames(datExpr0) <- toupper(colnames(datExpr0))
read.csv("blanca_hub_data_20minus.csv") -> fuzzy


## For interactomes that did NOT undergo fuzzy clustering (originally with size above 4 and below 21) ##

# Calculating individual correlation coefficients for each interaction
list_out <- NULL
for (j in unique(fuzzy$hub_id)) {
  library(WGCNA)
  print(j)
  subset(fuzzy, fuzzy$hub_id == j) -> test
  genes <- as.character(test$target)
  as.character(unique(test$source)) -> seed
  c(genes, seed) -> full
  sample_size <- round(nrow(datExpr0)*0.7, 0)
  sample(1:nrow(datExpr0), sample_size) -> rand_sample
  as.matrix(t(datExpr0[rand_sample,full])) -> x
  adjacency(t(x), type = "unsigned",
            power =1) -> df # now calculated outside
  list_apply <- sapply(genes, function(x) c("corr" = as.numeric(df[seed, x]),
                                            "target"= x, 
                                            "inter_id"= j, 
                                            "source"= seed, 
                                            "inter_size" = length(full)))
  as.data.frame(t(list_apply)) -> list_apply
  rbind(list_out, list_apply) -> list_out
}

# Calculating the average correlation coefficient for each interactome
df_means <- matrix(nrow = 0, ncol = 1)
for (j in unique(list_out$inter_id)) {
  print(j)
  list_out %>%
    filter(inter_id == j) -> module_list
  mean <- mean(as.numeric(data.matrix(module_list$corr)))
  for (i in 1:nrow(module_list)){
    df_means <- rbind(df_means, mean)
  }
  
}
list_out_means <- cbind(list_out, df_means)

# Save
write.csv(list_out, paste0(outdir, "truth_fuzzy_with_mod_under20.csv")) # individual interaction correlation coefficients
write.csv(list_out_means, paste0(outdir, "truth_fuzzy_with_mod_under20_MEANS.csv")) # interactome average correlation coefficient



## For interactomes that DID undergo fuzzy clustering (originally with size 21 and over) ##

read.csv("Fuzzy/gather_fanny.csv") -> fuzzy
fuzzy$interaction <- interaction(fuzzy$seed, fuzzy$source)

# Calculating individual correlation coefficients for each interaction
list_out <- NULL
for (j in unique(fuzzy$interaction)) {
  library(WGCNA)
  print(j)
  subset(fuzzy, fuzzy$interaction == j) -> test
  genes <- as.character(test$target)
  as.character(unique(test$seed)) -> seed
  c(genes, seed) -> full
  sample_size <- round(nrow(datExpr0)*0.7, 0)
  sample(1:nrow(datExpr0), sample_size) -> rand_sample
  as.matrix(t(datExpr0[rand_sample,full])) -> x
  adjacency(t(x), type = "unsigned",
            power =1) -> df # now calculated outside
  list_apply <- sapply(genes, function(x) c("corr" = as.numeric(df[seed, x]),
                                            "target"= x, 
                                            "inter_id"= j, 
                                            "source"= seed, 
                                            "inter_size" = length(full)))
  as.data.frame(t(list_apply)) -> list_apply
  rbind(list_out, list_apply) -> list_out
}

# Calculating the average correlation coefficient for each interactome
df_means <- matrix(nrow = 0, ncol = 1)
for (j in unique(list_out$inter_id)) {
  print(j)
  list_out %>%
    filter(inter_id == j) -> module_list
  mean <- mean(as.numeric(data.matrix(module_list$corr)))
  for (i in 1:nrow(module_list)){
    df_means <- rbind(df_means, mean)
  }
  
}
list_out_means <- cbind(list_out, df_means)

# Save
write.csv(list_out, paste0(outdir, "truth_fuzzy_with_mod_over20.csv")) # individual interaction correlation coefficients
write.csv(list_out_means, paste0(outdir, "truth_fuzzy_with_mod_over20_MEANS.csv")) # interactome average correlation coefficient


