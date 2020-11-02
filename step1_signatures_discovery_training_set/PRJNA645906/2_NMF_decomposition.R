# load required libraries and scripts
library("NMF")
library("nnls")
source("R/signatures.analysis.R")

# load data
load(file="results/contexts_matrix.RData")

# perform inference
signatures_decomposition = signatures.decomposition(x=contexts_matrix,K=1:6,nmf_runs=100,num_processes=10,seed=12345)

# save results
save(signatures_decomposition,file="results/signatures_decomposition.RData")
