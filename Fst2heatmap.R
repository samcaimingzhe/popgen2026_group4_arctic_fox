setwd("~/Desktop/AF/")
library(pheatmap)
library(tidyverse)
fst_df = read.csv('fst_values.csv')
pairs_list = strsplit(fst_df$X, "\\.vs\\.")
pops = sort(unique(unlist(pairs_list)))
n = length(pops)
fst_matrix = matrix(0, n, n, dimnames = list(pops, pops))
for(i in seq_along(pairs_list)) {
  pair <- pairs_list[[i]]
  fst_matrix[pair[1], pair[2]] = fst_df$fst_values[i]
  fst_matrix[pair[2], pair[1]] = fst_df$fst_values[i]
}
fst_matrix[fst_matrix<0]=0
pheatmap(fst_matrix, 
         display_numbers = TRUE,
         number_format = "%.3f",
         number_color = "black",
         fontsize_number = 10,
         cluster_rows = T, 
         cluster_cols = T,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Pairwise Fst between Populations",
         color = colorRampPalette(c("white", "blue"))(100))

