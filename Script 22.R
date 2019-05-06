
####### SCRIPT 22 #######

library("biomaRt")
library("dplyr")

## Load CSV data
whole_matrix <- read.csv("whole_matrix_changednames.csv", row.names = 1)
whole_matrix$X <- NULL
gene_names <- toupper(row.names(whole_matrix))


## Sample from the all-gene corr matrix and pass it through the API
random_genes <- matrix(nrow = 0, ncol = 2)
colnames(random_genes) <- c("gene_1","gene_2")
for (i in (1:549)){
  genes <- sample(gene_names, 2) # sample from 19,864 gene names
  random_genes <- rbind(random_genes, genes)
  genes <- paste(genes, collapse = "%0d")
  api_url <- paste("https://string-db.org/api/json/ppi_enrichment?identifiers=", genes, "&species=9606&required_score=400", sep="")
  response <- GET(api_url)
  response_r <- jsonlite::fromJSON(content(response, "text", encoding = "UTF8"), flatten = TRUE)
  print(response_r)
  response_list[[i]] <- response_r
}

random_api <- rbindlist(response_list, use.names = TRUE, idcol = T, fill = T)
results_total_df <- left_join(sig_interactions, results_truth_sig_api, by = c("nrow" = ".id"))

# manually add the results into a spreadsheet, number of known in one column, and number of sign in a different one
