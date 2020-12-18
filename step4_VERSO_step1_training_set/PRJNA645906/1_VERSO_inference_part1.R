# load required libraries and sources
library("ape")
library("parallel")
library("Rfast")
source("R/frontend.R")
source("R/inference.R")
source("R/utils.R")

# load variants data
load("results/clonal_variants.RData")

# select only frequent clonal variants
clonal_variants = clonal_variants[,names(which((colSums(clonal_variants,na.rm=TRUE)/nrow(clonal_variants))>0.03))]

# define the parameters of the grid search
alpha = rep(0.001,30)
beta = rep(0.001,30)

# perform inference
set.seed(12345)
inference_step1 = VERSO(D=clonal_variants,alpha=alpha,beta=beta,check_indistinguishable=TRUE,num_rs=1,num_iter=500000,n_try_bs=100000,num_processes=30,verbose=TRUE)

# save the results
save(inference_step1,file="results/inference_step1.RData")
