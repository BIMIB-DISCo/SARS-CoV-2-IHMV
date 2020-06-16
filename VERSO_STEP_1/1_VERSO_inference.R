# load required libraries and scripts
library("parallel")
library("Rfast")
source("R/frontend.R")
source("R/inference.R")
source("R/utils.R")

# load data
load("clonal_variants.RData")
D = list()
D[["Exp"]] = clonal_variants

# define parameters for grid search
alpha = NULL
beta = NULL
for(a in c(0.001,0.005,0.010,0.025,0.050,0.100,0.150,0.200)) {
    for(b in c(0.001,0.005,0.010,0.025,0.050,0.100,0.150,0.200)) {
        alpha = c(alpha,list(c(a)))
        beta = c(beta,list(c(b)))
    }
}

# perform inference
inference = VERSO(D=D,lik_w=1,alpha=alpha,beta=beta,keep_equivalent=TRUE,check_indistinguishable=TRUE,num_rs=500,num_iter=100000,n_try_bs=5000,marginalize=FALSE,num_processes=64,seed=12345,verbose=TRUE)

# save the results
save(inference,file="inference.RData")
