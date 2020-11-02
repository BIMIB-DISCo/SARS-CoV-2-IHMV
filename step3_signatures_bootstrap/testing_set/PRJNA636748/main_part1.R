# load required libraries and sources
library("glmnet")
library("nnls")
library("lsa")
library("parallel")
source("R/estimate.signatures.significance.R")

# read the results
load("RData/contexts_matrix.RData")
load("RData/signatures_decomposition.RData")

# perform a robust estimation of alpha coefficients
signatures_significance = signaturesSignificance(x=contexts_matrix,beta=signatures_decomposition$beta,cosine_thr=0.95,min_contribution=0.05,pvalue_thr=0.05,sparsify=FALSE,nboot=1000,num_processes=5,seed=44444,verbose=TRUE)

# save results
save(signatures_significance,file="results/signatures_significance.RData")
