
######### SCRIPT 21 #######

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("preprocessCore")

source("https://bioconductor.org/biocLite.R")
biocLite("GO.db")
biocLite("impute")
library(httr)
library(jsonlite)
library(dplyr)
library("data.table")
library("rlist")

# Load interactions with a significant mean truth outdegree
sig_interactions <- read.csv("significant_truth_matrix_39.csv", header = T)
sig_interactions$nrow <-  seq.int(nrow(sig_interactions))
response_list <- list()

# Loop querying the STRING database #
for (i in (1:549)){
  row <- sig_interactions[i,c(3,5)]
  genes <- unique(c(as.character(row$target), as.character(row$source)))
  genes <- paste(genes, collapse = "%0d")
  api_url <- paste("https://string-db.org/api/json/ppi_enrichment?identifiers=", genes, "&species=9606&required_score=400", sep="")
  response <- GET(api_url)
  response_r <- jsonlite::fromJSON(content(response, "text", encoding = "UTF8"), flatten = TRUE)
  print(response_r)
  response_list[[i]] <- response_r
}

results_truth_sig_api <- rbindlist(response_list, use.names = TRUE, idcol = T)
results_total_df <- left_join(sig_interactions, results_truth_sig_api, by = c("nrow" = ".id"))
# 10/549 have PPI enrichment p-values below 0.05 - 1.8%

write.csv(results_total_df, "STRING_results_39.csv", row.names = FALSE)

