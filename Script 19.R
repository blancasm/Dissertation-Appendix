
###### SCRIPT 19 - fuzzy clustering of hubs identified from first outdegree with size > 20 ######

### nb. Script by Dr. Tina Baker, but I edited parts of it ###

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("preprocessCore")

source("https://bioconductor.org/biocLite.R")
biocLite("GO.db")
biocLite("impute")

library(WGCNA)
library(igraph)
library(randomcoloR)
library(reshape2)
library(rsgcc)
library(ggpubr)
library(ggplot2)
library(networkD3)
library(factoextra)
library(tidyverse)
library(nnet)
library(ape)
library(dendextend)
library(gplots) 
library(tidyr)
library(cluster)
library(caret)

source("Fuzzy/flattenCorrMatrix.R") ## TO CHANGE  ##
options(stringsAsFactors = FALSE)
################ READ IN RNASEQ DATA ############################################
# Read in the data
file = "Fuzzy/tpm_0.001_log2_merge_dedup_symb_ACT_LTB_PPD_for_eqtl.csv" 
datExpr0 = read.csv(file = file, row.names = "X")
outdir <- "Fuzzy/" ## TO CHANGE  ##
################ READ IN EDGE FILE #############################################
# read in edge file
read.csv("blanca_hub_data_20plus.csv") -> edges 
#which(grepl(pattern = '[a-z]', x = edges$source))
toupper(edges$source) -> edges$source
toupper(edges$target) -> edges$target
# remove self loops 
caret::nearZeroVar(t(edges[, c("source", "target")])) -> select_rows
edges[-select_rows,] -> edges
length(unique(edges$source))
# Original Edge table connections
as.data.frame(table(edges$source)) -> total_inter
sum(total_inter$Freq) 
base::subset(edges, edges$source %in% colnames(datExpr0) & edges$target %in% colnames(datExpr0)) -> edges_match
# Matched edge table connections
as.data.frame(table(edges_match$source)) -> total_match
sum(total_match$Freq) 
################## SUBSET TO ONLY 21 + SIZE MODULES ############################
base::subset(total_match, total_match$Freq >=21) %>%
  droplevels() -> total_match_inter
rownames(total_match_inter) <- NULL

which(grepl(colnames(datExpr0), pattern = "[A-Z]"))
which(grepl(edges_match$target, pattern = "[a-z]"))

subset(edges_match, edges_match$source == "CPM") -> cpm
length(unique(cpm$target))
#################### START LOOP ################################################
fuzzy_clust_min20 <- function(total_match_inter, edges_match, outdir){
  gc()
  df_fanny_gather_out <- NULL
  results_fuzzy_out <- NULL
  list_out <- NULL
  mean_out = NULL
  for (j in 1:length(total_match_inter$Var1)) {
    total_match_inter$Var1[j] -> i
    total_match_inter$Freq[j] -> deg
    subset(edges_match, edges_match$source == i) -> test
    genes <- unique(test$target)
    as.character(i) -> seed ; c(genes, seed) -> full
    print(seed)
    sample_size <- round(nrow(datExpr0)*0.7, 0)
    sample(1:nrow(datExpr0), sample_size) -> rand_sample
    ############## Select random sample and scale #####
    as.matrix(t(datExpr0[rand_sample,])) -> datExpr0_sample
    scaledata <- t(scale(t(datExpr0))) # Centers and scales data.
    scaledata <- scaledata[complete.cases(scaledata),]
    ############### with seed ####################
    as.matrix(t(scaledata[rand_sample,full])) -> x 
    adjacency(t(x), type = "unsigned", # unsigned or unsigned
              power =1) -> df
    flattenCorrMatrix(df) -> pcc
    mean(abs(as.numeric(pcc$cor)), na.rm = TRUE) -> avg
    ############### without seed ##################
    as.matrix(t(scaledata[rand_sample,genes])) -> x2 ### change here
    adjacency(t(x2), type = "unsigned", # unsigned or unsigned
              power =1) -> df2 
    flattenCorrMatrix(df2) -> pcc2
    mean(abs(as.numeric(pcc2$cor)), na.rm = TRUE) -> avg2
    ############### output summary ###############
    mean <- data.frame("sample_num" = i, "mod_size" = length(full), "with_seed" = avg, 
                       "without_seed" = avg2, seed)
    rbind(mean_out, mean) -> mean_out
    ################# fuzzy clustering ###########
    # n_clust <- fviz_nbclust(x2, kmeans,  #fuzzy k means (done on distance matrix)
    #                         method = "silhouette", 
    #                         k.max = nrow(x2)-1) # need kmax here -1 otherwise error ?? 
    # k <- which.is.max(n_clust$data$y); print(paste0("k = ",k)) # setting k
    ################ fuzzy ###################
    d <- dist(df2) # only take without seed
    fannyy <- fanny(d, k=round((nrow(x2)/10), 0), 
                    memb.exp = 1.2,
                    maxit = 1000)
    round(fannyy$membership, 2); fannyyMA <- round(fannyy$membership, 2) > 0.10 
    print(apply(fannyyMA, 1, function(x) paste(which(x), collapse="_")))
    fannyyG <- fannyy$clustering 
    ################# format membership #########
    as.data.frame(fannyyMA) %>%  
      rownames_to_column(., 'target') %>% 
      gather(., "source", "value", -"target") %>%
      base::subset(., value == "TRUE") %>%
      mutate("seed" = seed)-> df_fanny_gather
    as.data.frame(fannyyG) %>%  ## add max group to force no group
      rownames_to_column(., 'target') %>%
      mutate("seed" = seed) %>%
      mutate("source" = paste0("V",fannyyG)) %>%
      mutate("value" = "TRUE")->  group_force
    group_force$fannyyG <- NULL
    rbind(group_force, df_fanny_gather) -> df_fanny_gather_force
    duplicated(df_fanny_gather_force)
    unique(df_fanny_gather_force) -> df_fanny_gather_force; 
    rownames(df_fanny_gather_force) <- NULL
    as.data.frame(fannyyG) %>%  
      rownames_to_column(., 'target') -> group
    melt(x2) -> melted
    left_join(melted, group, by = c("Var1" = "target")) -> joined
    ############### plot groups ###################
    k_plot <- ggplot(data = joined, 
                     aes(x = Var2, y = value, group = Var1)) +#colour = as.factor(fannyyG)
      geom_line() +
      facet_wrap(as.factor(fannyyG)) +
      theme(axis.text.x = element_blank())
    ppi <- 300
    png(paste0(outdir,seed,"_ks.png"), width=8*ppi, height=8*ppi, res=ppi, pointsize =14)
    plot(k_plot)
    dev.off() -> off
    ############## make full results file  #########
    results_out <- NULL
    df_fanny_gather_force -> df_fanny_gather
    for (z in unique(df_fanny_gather$source)){
      mod <- base::subset(df_fanny_gather, df_fanny_gather$source == z)
      mod_genes <- as.character(mod$target)
      as.matrix(t(scaledata[rand_sample ,c(mod_genes,seed)])) -> mod_df
      as.matrix(t(scaledata[rand_sample, mod_genes])) -> mod_df2
      adjacency(t(mod_df), type = "unsigned", # unsigned or unsigned
                power =1) -> mod_df_corr
      adjacency(t(mod_df2), type = "unsigned", # unsigned or unsigned
                power =1) -> mod_df_corr2
      flattenCorrMatrix(mod_df_corr) -> mod_pcc
      mean(abs(as.numeric(mod_pcc$cor)), na.rm = TRUE) -> avg
      flattenCorrMatrix(mod_df_corr2) -> mod_pcc2
      mean(abs(as.numeric(mod_pcc2$cor)), na.rm = TRUE) -> avg2
      results <- data.frame("clust_num" = z, "clust_size" = length(mod_genes), "with_seed" = avg,
                            "without_seed" = avg2, seed)
      results_out <- rbind(results_out, results)
    }
    results_fuzzy_out <- rbind(results_fuzzy_out, results_out)
    df_fanny_gather_out <- rbind(df_fanny_gather_out, df_fanny_gather)
    ################## plot fuzzy clust #########
    #plot_fuzzy_clust_FUN(df_fanny_gather = df_fanny_gather, outdir = outdir, seed = seed)
    n_clust <- NULL
    fannyy <- NULL
    ################## Rbind fanny out #################
    list_out <- rbind(list_out, df_fanny_gather)
  }# end of main for loop
  write.csv(df_fanny_gather_out, paste0(outdir, "gather_fanny.csv"), row.names = FALSE)
  write.csv(results_fuzzy_out, paste0(outdir, "results_fuzzy_out.csv"), row.names = FALSE)
  write.csv(mean_out, paste0(outdir, "mean_out.csv"), row.names = FALSE)
  return(results_fuzzy_out)
}# end of function


fuzzy_clust_min20(total_match_inter= total_match_inter, 
                  edges_match = edges_match,  
                  outdir= outdir) -> results_fuzzy_out


