# Our Population Genetics Group Project: Arctic Fox
Final report is here (PDF format)
[Report Link](https://github.com/samcaimingzhe/popgen2026_group4_arctic_fox/blob/main/Project_report_Group4.pdf)
## Pre-processing
> All codes are avalible in this repository.
> Datasets are avalible in my path /home/popgenmsc26user21/github, you may cp them if you need.

Copy the data to your path.
```
cp -r /course/popgenmsc26/projects/arctic_fox/data/* ./
```
Change the name of chromosomes into integer with python script.
```python
with open("AF.imputed.thin.0.bim", "r") as infile:
    content = infile.read()

for i in range(1, 25):
    acc = f"NC_{54823+i:06d}.1"
    content = content.replace(acc, str(i))

with open("AF.imputed.thin.bim", "w") as outfile:
    outfile.write(content)

print("Done!")
```
```
cp AF.imputed.thin.bim AF.imputed.thin.0.bim
python3 chr.py
```
## PCA
Obtain PCA results by plink. 10 PCs will be calculted though only first 2 are useful.
```
plink --bfile AF.imputed.thin --freq --chr-set 24 --out AF_freqs
plink --bfile AF.imputed.thin --read-freq AF_freqs.frq --pca 10 --chr-set 24 --out AF_pca
```
Example code to download the pca data files to your local for plotting.

NB: remember to change the user number!
```
mkdir -p ~/Desktop/AF
scp "popgenmsc26user00@emily.popgen.dk:~/github/AF_pca.*" ~/Desktop/AF
scp "popgenmsc26user00@emily.popgen.dk:~/github/sample_popinfo.tsv" ~/Desktop/AF
```
```R
setwd('~/Desktop/AF') # make a AF folder on your desktop
data = read.table('AF_pca.eigenvec')
popinfo = read.table('sample_popinfo.tsv',header = T)
colnames(data) <- c('fid', 'iid', 
                    'pc1', 'pc2', 'pc3','pc4', 'pc5',
                    'pc6','pc7', 'pc8', 'pc9','pc10')
data$fid = popinfo$Region

library(plotly)
plot_ly(data, 
        x = ~pc1, 
        y = ~pc2, 
        z = ~pc3, 
        color = ~fid,     
        colors = "Set1",  
        type = 'scatter3d', 
        mode = 'markers',
        marker = list(size = 4, opacity = 0.8),
        text = ~iid)
```
<img width="729" height="476" alt="pca" src="https://github.com/user-attachments/assets/82f2f873-3922-41ac-bc99-f4803806e810" />
So we observed that ArcFoxSamp_128, ArcFoxSamp_144, ArcFoxSamp_152 may be the 3 most recent immigrants.

## LD Pruning
Admixture assume the SNPs are independent to each others, so we have to do LD prunning.
```
# remove SNPs with missing genotypes
plink --bfile AF.imputed.thin --geno 0 --make-bed --out AF.complete
# remove low MAF SNPs
plink --bfile AF.complete --maf 0.05 --make-bed --out AF.frq_normal
# LD pruning
plink --bfile AF.frq_normal --indep-pairwise 50 10 0.1 --out AF.LD
plink --bfile AF.frq_normal --extract AF.LD.prune.in  --make-bed --out AF.clean
```
Then we can check how many SNPs are filtered out during each step.
```
wc -l *.bim|sort -r
```
```
 2906282 total
  816284 AF.imputed.thin.bim
  816284 AF.imputed.thin.0.bim
  719391 AF.complete.bim
  498777 AF.frq_normal.bim
   55546 AF.clean.bim
```
## Admixture \& Evaladmix processing
Then we can run the admixture by K loop and run loop, try to run it in `screen`:
```
screen -S admix
```
Inside admix, we run this double loop commands:
```
for K in {2..7}; do
  for i in {1..100}; do
    echo "running K=${K} & run=${i}"
    admixture -s ${i} AF.clean.bed ${K} > AF.clean_K${K}_run${i}.log
    cp AF.clean.${K}.Q AF.clean_K${K}_run${i}.Q
    cp AF.clean.${K}.P AF.clean_K${K}_run${i}.P          
  done
done
```
You may use `control + A + D` to exit the screen and log out the server, because admixture requires a while to run.
We can continue to run Fst test. Please download the files to local.
```
scp "popgenmsc26user00@emily.popgen.dk:~/github/AF.imputed.thin.???" ~/Desktop/AF
wget https://raw.githubusercontent.com/samcaimingzhe/popgen2026_group4_arctic_fox/main/modified_popinfo.tsv -P ~/Desktop/AF
```
## $F_{st}$
Here is the Rscript for $F_{st}$ calculation and heatmap visualization, and if you don't want to calculate locally, we also provide the results data. 

Here is the result of $F_{st}$ calculation:
```
wget https://raw.githubusercontent.com/samcaimingzhe/popgen2026_group4_arctic_fox/refs/heads/main/fst_values.csv
```
Here is the modified popinfo:
```
wget https://raw.githubusercontent.com/samcaimingzhe/popgen2026_group4_arctic_fox/main/modified_popinfo.tsv
```
The $F_{st}$ calculation:
```R
setwd("~/Desktop/AF/")
WC84<-function(x,pop){
  n<-table(pop)
  npop<-nrow(n)
  n_avg<-mean(n)
  N<-length(pop)
  p<-apply(x,2,function(x,pop){tapply(x,pop,mean)/2},pop=pop)
  p_avg<-as.vector(n%*%p/N )
  s2<-1/(npop-1)*(apply(p,1,function(x){((x-p_avg)^2)})%*%n)/n_avg
  h_avg<-apply(x==1,2,sum)/N
  n_c<-1/(npop-1)*(N-sum(n^2)/N)
  a <-n_avg/n_c*(s2-(p_avg*(1-p_avg)-(npop-1)*s2/npop-h_avg/4)/(n_avg-1))
  b <- n_avg/(n_avg-1)*(p_avg*(1-p_avg)-(npop-1)*s2/npop-(2*n_avg-1)*h_avg/(4*n_avg))
  c <- h_avg/2
  F <- 1-c/(a+b+c)
  theta <- a/(a+b+c)
  f <- 1-c(b+c)
  theta_w<-sum(a)/sum(a+b+c)
  list(F=F,theta=theta,f=f,theta_w=theta_w,a=a,b=b,c=c,total=c+b+a)
}

library(snpStats)
data <- read.plink("AF.imputed.thin")
geno <- matrix(as.integer(data$genotypes),nrow=nrow(data$genotypes))
geno[geno==0] <- NA
geno <- geno - 1
g <- geno[,complete.cases(t(geno))]
dim(geno)
dim(g)

popinfo <- read.table("modified_popinfo.tsv", stringsAsFactors=F, header = T)
region <- unique(popinfo$Region)
sapply(region, function(x) popinfo$Sample[popinfo$Region == x])
region_pairs <- t(combn(region, 2))

fsts <- apply(region_pairs, 1, function(x) WC84(g[popinfo$Region %in% x,], 
                                                popinfo$Region[popinfo$Region %in% x]))

names(fsts) <- apply(region_pairs, 1, paste, collapse=".vs.")
lapply(fsts, function(x) x$theta_w)
fst_values <- sapply(fsts, function(x) x$theta_w)
fst_values.df = data.frame(fst_values)
write.csv(fst_values.df,'fst_values.csv')
```
The heatmap visualization:
```R
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
```
You may get a heatmap shows the pairwise $F_{st}$ between populations with cluster:
<img width="1114" height="1012" alt="Fst" src="https://github.com/user-attachments/assets/b3f0c5b2-9d79-4af0-9cbe-f928975d79af" />

## Admixture \& Evaladmix analysis
After a long time we may complete the **Admixure** (appoximately 9 hours):
```
for K in {2..7}; do
    for i in {1..100}; do
        log=$(grep '^Loglikelihood' AF.clean_K${K}_run${i}.log | awk '{print $2}')
        echo -e "${K}\t${i}\t${log}\tAF.clean_K${K}_run${i}.log" >> Loglikelihoods.log
    done
done

sort -t$'\t' -k1,1 -k3,3gr Loglikelihoods.log | awk -F'\t' 'count[$1]++ < 3' > top3_log_likelihood.log
sort -t$'\t' -k1,1 -k3,3gr Loglikelihoods.log | awk -F'\t' 'count[$1]++ < 1' > top1_log_likelihood.log
```
You may see K=2-5 are converged, but all top 1 are going to be **Evaladmix**:
```
cat top3_log_likelihood.log

2	27	-2820433.464874	AF.clean_K2_run27.log
2	5	-2820433.464874	AF.clean_K2_run5.log
2	50	-2820433.464874	AF.clean_K2_run50.log
3	68	-2755472.370008	AF.clean_K3_run68.log
3	1	-2755472.370010	AF.clean_K3_run1.log
3	35	-2755472.370010	AF.clean_K3_run35.log
4	93	-2712268.360483	AF.clean_K4_run93.log
4	30	-2712268.360519	AF.clean_K4_run30.log
4	21	-2712268.360523	AF.clean_K4_run21.log
5	63	-2672566.870710	AF.clean_K5_run63.log
5	13	-2672566.871107	AF.clean_K5_run13.log
5	93	-2672566.871704	AF.clean_K5_run93.log
6	25	-2633122.089118	AF.clean_K6_run25.log
6	67	-2633122.113995	AF.clean_K6_run67.log
6	48	-2633457.568546	AF.clean_K6_run48.log
7	78	-2594992.707457	AF.clean_K7_run78.log
7	23	-2595165.099862	AF.clean_K7_run23.log
7	63	-2595276.661517	AF.clean_K7_run63.log
```
By selecting:
```
cat top1_log_likelihood.log

2	27	-2820433.464874	AF.clean_K2_run27.log
3	68	-2755472.370008	AF.clean_K3_run68.log
4	93	-2712268.360483	AF.clean_K4_run93.log
5	63	-2672566.870710	AF.clean_K5_run63.log
6	25	-2633122.089118	AF.clean_K6_run25.log
7	78	-2594992.707457	AF.clean_K7_run78.log
```
Copy the P file to bestAdmix for further download and plots, also run Evaladmix at the same time:
```
mkdir -p bestAdmix
cp /course/popgenmsc26/exercises/structure/evalAdmix ./bestAdmix
cp /course/popgenmsc26/exercises/structure/visFuns.R ./bestAdmix

while read -r line; do
    K=$(echo "$line" | awk '{print $1}')
    file_raw=$(echo "$line" | cut -f4)
    file_prefix=${file_raw%.log}
    cp "${file_prefix}.P" ./bestAdmix/
    cp "${file_prefix}.Q" ./bestAdmix/
    echo "Copied ${file_prefix} to ./bestAdmix"
    ./evalAdmix \
        -plink AF.clean \
        -fname bestAdmix/${file_prefix}.P \
        -qname bestAdmix/${file_prefix}.Q \
        -o ./bestAdmix/K${K}.corres.txt

    echo "Finished evalAdmix for K=$K. Result saved to ./bestAdmix/K${K}.corres.txt"
done < top1_log_likelihood.log
```
By downloading the Q files to local:
```
scp "popgenmsc26user00@emily.popgen.dk:~/github/bestAdmix/*.Q" ~/Desktop/AF
scp "popgenmsc26user00@emily.popgen.dk:~/github/bestAdmix/*corres.txt" ~/Desktop/AF
scp "popgenmsc26user00@emily.popgen.dk:~/github/bestAdmix/visFuns.R" ~/Desktop/AF
```

Plot the admixture and evaladmix results:
```R
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

plot_evaladmix('K2.corres.txt', 2)
plot_evaladmix('K3.corres.txt', 3)
plot_evaladmix('K4.corres.txt', 4)
plot_evaladmix('K5.corres.txt', 5)
plot_evaladmix('K6.corres.txt', 6)
plot_evaladmix('K7.corres.txt', 7)
```
This is the result of **Evaladmix**, we have the K=5 as the best K for our data, but from a pratical perspective, Taymr does not neccessarily to be 0, so K=4 maybe also enough? (Just a guess, we didn't ask Xiaodong yet).
<img width="3121" height="1677" alt="evaladmix" src="https://github.com/user-attachments/assets/655a8a8f-7e6b-4399-adab-c13e5abe714b" />
Here is the result of **Admixture**:
[admixture_plot.pdf](https://github.com/user-attachments/files/26102004/admixture_plot.pdf)


And when we focus on 3 immigrants when K=5, these patterns are likely to tell us where did they come from, but we don't have the data to inference this, which is we would like to go further:
<img width="895" height="340" alt="截屏2026-03-19 上午12 02 09" src="https://github.com/user-attachments/assets/815d2bc2-e3ec-49da-99fc-e85df22e1927" />

## Genetic Diversity
```
mkdir -p subpop
wget https://raw.githubusercontent.com/samcaimingzhe/popgen2026_group4_arctic_fox/refs/heads/main/freq_popinfo.tsv
for p in Qanisartuut Zackenberg Scoresbysund Kangerlussuaq Scoresbysund_immigrant Canada_Siberia
do
    echo "Processing $p..."
    awk -v reg="$p" '$1==reg{print"0",$2}' freq_popinfo.tsv > keep_${p}.txt

    plink --bfile AF.imputed.thin \
          --keep keep_${p}.txt \
          --mac 1 \
          --geno 0.5 \
          --make-bed \
          --allow-no-sex \
          --out subpop/$p \
          --silent
          
    plink --bfile subpop/$p --freq --out subpop/$p
    plink --bfile subpop/$p --freq --out subpop/$p
done
```
Download these `.frq` files to local:
```
scp "popgenmsc26user21@emily.popgen.dk:~/github/subpop/*.frq" ~/Desktop/AF
```
Draw the boxplot:
```
setwd("~/Desktop/AF/")
library(ggplot2)
library(dplyr)
files <- c("Canada_Siberia.frq", "Kangerlussuaq.frq", "Qanisartuut.frq", 
           "Scoresbysund.frq", "Scoresbysund_immigrant.frq", "Zackenberg.frq")

data_het = function(FILE){
  data = read.table(FILE, header = T)
  data$het = 2 * data$MAF * (1 - data$MAF)

  filename = basename(FILE)
  name = sub("\\.frq$", "", filename)

  result = data.frame(
    name = name,
    het = data$het
  )
  
  return(result)
}

all_data = rbind(data_het("Canada_Siberia.frq"),
data_het("Kangerlussuaq.frq"),
data_het("Qanisartuut.frq"),
data_het("Scoresbysund.frq"),
data_het("Scoresbysund_immigrant.frq"),
data_het("Zackenberg.frq")
)

ggplot(all_data, aes(x = name, y = het, fill = name, group = name)) +
  coord_flip() +
  geom_boxplot() +
  labs(y = 'Heterozygosity', x = '') +
  theme_light() +
  theme(legend.position = 'none')
```
<img width="790" height="426" alt="het" src="https://github.com/user-attachments/assets/e81b66e2-ddff-45fd-bc26-60d417fdd2ec" />

## $\pi$, $\theta_w$, and Tajima's D
```
setwd("~/Desktop/AF/")

# Function for expected heterozygosity
het <- function(x) { 2 * x * (1 - x) }

# Load one population (example: Kangerlussuaq)
load_pop <- function(prefix) {
  frq <- read.table(paste0(prefix, ".frq"), h = TRUE)
  bim <- read.table(paste0(prefix, ".bim"), h = FALSE)
  
  frq$het <- het(frq$MAF)
  frq$pos <- bim$V4
  frq$CHR <- frq$CHR
  
  return(frq)
}

# Estimate sample size
getN <- function(x) {
  round(median(x$NCHROBS, na.rm = TRUE) / 2)
}

# Per chromosome stats
calc_chr_stats <- function(df) {
  n <- getN(df)
  S <- nrow(df)
  
  if (n < 2 || S == 0) return(NULL)
  
  a1 <- sum(1 / (1:(n - 1)))
  a2 <- sum(1 / ((1:(n - 1))^2))
  
  b1 <- (n + 1) / (3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))
  
  c1 <- b1 - 1 / a1
  c2 <- b2 - (n + 2) / (a1 * n) + a2 / (a1^2)
  
  e1 <- c1 / a1
  e2 <- c2 / (a1^2 + a2)
  
  span <- max(df$pos) - min(df$pos)
  
  if (span == 0) return(NULL)
  
  k <- sum(df$het, na.rm = TRUE)
  
  pi_chr <- k / span
  theta_chr <- S / (a1 * span)
  tajima_chr <- (k - S / a1) / sqrt(e1 * S + e2 * S * (S - 1))
  
  return(data.frame(
    CHR = unique(df$CHR),
    span = span,
    pi = pi_chr,
    thetaW = theta_chr,
    TajimaD = tajima_chr
  ))
}

# Combine across chromosomes
combine_stats <- function(df) {
  # weighted mean for pi and thetaW
  pi_total <- sum(df$pi * df$span) / sum(df$span)
  theta_total <- sum(df$thetaW * df$span) / sum(df$span)
  
  # Tajima's D mean and SD
  tajima_mean <- mean(df$TajimaD, na.rm = TRUE)
  tajima_sd <- sd(df$TajimaD, na.rm = TRUE)
  
  return(list(
    pi = pi_total,
    thetaW = theta_total,
    TajimaD_mean = tajima_mean,
    TajimaD_sd = tajima_sd
  ))
}

# Run for one population
analyze_pop <- function(prefix) {
  df <- load_pop(prefix)
  
  chr_list <- split(df, df$CHR)
  
  chr_stats <- do.call(rbind, lapply(chr_list, calc_chr_stats))
  
  combined <- combine_stats(chr_stats)
  
  return(list(
    chr_stats = chr_stats,
    combined = combined
  ))
}

# secondary level Populations
pop_names <- c("Kangerlussuaq", "Qanisartuut", 
               "Scoresbysund", "Scoresbysund_immigrant",
               "Zackenberg", "Canada_Siberia")

results <- lapply(pop_names, analyze_pop)
names(results) <- pop_names

# Extract final values
all_pi <- sapply(results, function(x) x$combined$pi)
all_thetaW <- sapply(results, function(x) x$combined$thetaW)
all_tajimaD <- sapply(results, function(x) x$combined$TajimaD_mean)

# Show results
all_pi
all_thetaW
all_tajimaD

# Plot pi
par(mar = c(7,10,4,5))
barplot(all_pi,
        names.arg = pop_names,
        las = 2,
        horiz = TRUE,
        xlim = c(0,0.00007),
        cex.names = 0.8,
        xlab = expression(pi))

# Plot thetaW
par(mar = c(7,10,4,2))
barplot(all_thetaW,
        names.arg = pop_names,
        las = 2,
        horiz = TRUE,
        xlim = c(0,2e-04),
        cex.names = 0.8,
        xlab = expression(theta[w]),
        mgp = c(4, 1, 0))

# Extract mean and SD
all_tajimaD <- sapply(results, function(x) x$combined$TajimaD_mean)
all_tajima_sd <- sapply(results, function(x) {
  sd(x$chr_stats$TajimaD, na.rm = TRUE)
})

pop_names <- names(results)

xlim_range <- c(min(all_tajimaD - all_tajima_sd, na.rm = TRUE), 0.5)
par(mar = c(5, 10, 4, 2))

bp <- barplot(all_tajimaD,
              names.arg = pop_names,
              las = 1,     
              horiz = TRUE,
              cex.names = 0.8,
              xlab = "Tajima's D",
              xlim = xlim_range)

arrows(x0 = all_tajimaD - all_tajima_sd,
       y0 = bp,                         
       x1 = all_tajimaD + all_tajima_sd,
       y1 = bp,                         
       angle = 90,
       code = 3,
       length = 0.05)

abline(v = 0, lty = 2, col = "red", lwd = 2)
```
<img width="948" height="1012" alt="pi" src="https://github.com/user-attachments/assets/a701d1af-e06c-49e8-83cc-26b995bc5503" />
<img width="948" height="1012" alt="thetaw" src="https://github.com/user-attachments/assets/d061e616-f145-4a8a-8129-82fa0f075814" />
<img width="946" height="878" alt="tajimaD" src="https://github.com/user-attachments/assets/8fc7c5b7-f089-4573-9df9-1bf22ee55377" />








### TO BE UPDATED ...
