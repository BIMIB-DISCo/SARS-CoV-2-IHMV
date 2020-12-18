# load required libraries and sources
library("ape")
library("parallel")
library("Rfast")
source("R/frontend.R")
source("R/inference.R")
source("R/utils.R")
source("R/VERSO_v2.R")
source("R/VERSOphylo.R")

# load variants data
load("results/clonal_variants.RData")

# select only frequent clonal variants
clonal_variants = clonal_variants[,names(which((colSums(clonal_variants,na.rm=TRUE)/nrow(clonal_variants))>0.03))]
load(file="results/inference_step1.RData")

# define the parameters of the grid search
alpha = NULL
beta = NULL
for(a in c(0.001,0.010,0.050)) {
    for(b in c(0.001,0.010,0.050)) {
        alpha = c(alpha,a)
        beta = c(beta,b)
    }
}
alpha = c(alpha,alpha,alpha)
beta = c(beta,beta,beta)

# set initialization
initialization = inference_step1$B
clonal_variants = clonal_variants[,-which(colnames(clonal_variants)=="28883_G_C")]
colnames(clonal_variants)[which(colnames(clonal_variants)=="28882_G_A")] = "28882_G_A|28883_G_C"
for(i in 2:ncol(initialization)) {
    colnames(initialization)[i] = as.character(which(colnames(clonal_variants)==colnames(initialization)[i]))
}
rownames(initialization)[2:nrow(initialization)] = 1:(nrow(initialization)-1)
rownames(initialization)[1] = "r"
colnames(initialization)[1] = "r"

# perform inference
set.seed(54321)
inference_step2 = VERSO_v2(D=clonal_variants,alpha=alpha,beta=beta,check_indistinguishable=TRUE,num_rs=1,num_iter=500000,n_try_bs=100000,num_processes=27,verbose=TRUE,initialization=initialization)

# save the results
save(inference_step2,file="results/inference_step2.RData")
