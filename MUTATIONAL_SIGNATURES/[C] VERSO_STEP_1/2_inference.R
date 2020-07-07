# load required libraries and scripts
library("parallel")
library("Rfast")
source("R/frontend.R")
source("R/inference.R")
source("R/utils.R")

# load data
load("raw_data/clonal_variants.RData")
D = list()
D[["Exp"]] = clonal_variants

# define parameters for grid search
alpha = NULL
beta = NULL
for(a in c(0.001,0.005,0.010,0.015)) {
    for(b in c(0.001,0.005,0.010,0.015)) {
        alpha = c(alpha,list(c(a)))
        beta = c(beta,list(c(b)))
    }
}

# perform inference
inference = VERSO(D=D,lik_w=1,alpha=alpha,beta=beta,keep_equivalent=TRUE,check_indistinguishable=TRUE,num_rs=100,num_iter=10000,n_try_bs=1000,marginalize=FALSE,num_processes=16,seed=3650499,verbose=TRUE)

# save the results
save(inference,file="results/inference.RData")
