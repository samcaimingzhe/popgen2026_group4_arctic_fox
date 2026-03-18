library(ggplot2)
library(tidyr)
library(dplyr)

plot_ancestry_ggplot = function(DATA_NAME, TITEL, K) {
  snp = read.table(DATA_NAME)
  colnames(snp) = paste0("V", 1:K)
  
  snp$Sample = popinfo$Sample
  snp$Region = popinfo$Region
  
  region_order = c("Qanisartuut", "Scoresbysund", "Zackenberg", "Kangerlussuaq", 
                   "Bylot_island", "Karrak_lake", "Scoresbysund_immigrant", "Taymyr")
  
  snp$Region = factor(snp$Region, levels = region_order)
  
  snp = snp %>% arrange(Region, Sample)
  snp$Sample = factor(snp$Sample, levels = unique(snp$Sample))
  
  df_long = snp %>%
    pivot_longer(cols = starts_with("V"), 
                 names_to = "Ancestor", 
                 values_to = "Proportion")
  
  colors = if(K <= 2) c("#E41A1C", "#377EB8")[1:K] else RColorBrewer::brewer.pal(min(K, 9), "Set1")[1:K]
  
  p = ggplot(df_long, aes(x = Sample, y = Proportion, fill = Ancestor)) +
    geom_col(width = 0.9) + 
    facet_grid(~Region, scales = "free_x", space = "free_x") + 
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = c(0, 0)) + 
    theme_minimal() +
    labs(title = TITEL, x = NULL, y = "Proportion") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5), 
      axis.ticks.x = element_line(size = 0.2),
      panel.spacing = unit(0.1, "lines"), 
      legend.position = "none",
      panel.grid = element_blank(),
      strip.text.x = element_blank(),
      strip.background = element_blank(),
    )
  
  return(p)
}

popinfo = read.table("modified_popinfo.tsv",header = T)
fam = read.table("AF.imputed.thin.fam",header = F)

library(patchwork)
p2 = plot_ancestry_ggplot("AF.clean_K2_run27.Q", "K=2", 2)
p3 = plot_ancestry_ggplot("AF.clean_K3_run68.Q", "K=3", 3)
p4 = plot_ancestry_ggplot("AF.clean_K4_run93.Q", "K=4", 4)
p5 = plot_ancestry_ggplot("AF.clean_K5_run63.Q", "K=5", 5)
p6 = plot_ancestry_ggplot("AF.clean_K6_run25.Q", "K=6", 6)
p7 = plot_ancestry_ggplot("AF.clean_K7_run78.Q", "K=7", 7)

final_plot = (p2 / p3 / p4 / p5 / p6 / p7) + 
  plot_annotation(title = "Admixture Analysis (K=2-7)",
                  theme = theme(plot.title = element_text(size = 20, hjust = 0.5)))

ggsave("admixture_plot.pdf", plot = final_plot, width = 7, height = 16)


source("visFuns.R")
library(ggplot2)
library(dplyr)
library(reshape2)

popinfo <- read.table("modified_popinfo.tsv",header=T)

pop = as.vector(popinfo[,1])
plot_evaladmix = function(FILE,K){
  r <- as.matrix(read.table(FILE))
  plotCorRes(cor_mat = r, pop = pop, 
             title = paste0("Correlation of residuals (K=", K,")"), 
             max_z=0.15, min_z=-0.15,
             cex.main=1, cex.lab=0.7, cex.legend=1,
             rotatelabpop=90,
             pop_labels = c(T,F),adjlab = 0.2)
  }

par(mfrow = c(2, 3))
plot_evaladmix('K2.corres.txt', 2)
plot_evaladmix('K3.corres.txt', 3)
plot_evaladmix('K4.corres.txt', 4)
plot_evaladmix('K5.corres.txt', 5)
plot_evaladmix('K6.corres.txt', 6)
plot_evaladmix('K7.corres.txt', 7)

