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
