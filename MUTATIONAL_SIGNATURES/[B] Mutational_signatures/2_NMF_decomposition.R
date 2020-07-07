# load required libraries and scripts
library("NMF")
library("nnls")
source("R/signatures.analysis.R")

# load data
load(file="results/contexts_matrix.RData")

# perform inference
signatures_decomposition = signatures.decomposition(x=contexts_matrix,K=1:6,nmf_runs=1000,num_processes=10,seed=12345)
signatures_cross_validation = signatures.cross.validation(x=(contexts_matrix/rowSums(contexts_matrix)),beta=signatures_decomposition$beta,cross_validation_entries=0.01,cross_validation_iterations=5,cross_validation_repetitions=1000,num_processes=10,seed=54321)
signatures_results = list(decomposition=signatures_decomposition,cross_validation=signatures_cross_validation)

# save results
save(signatures_results,file="results/signatures_results.RData")
