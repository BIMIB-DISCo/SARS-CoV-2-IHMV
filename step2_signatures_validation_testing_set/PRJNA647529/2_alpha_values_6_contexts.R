# load results
load("results/contexts_matrix.RData")
load("training/signatures_decomposition.RData")
observed_counts = contexts_matrix
beta = signatures_decomposition$beta[[3]]

# perform fit
set.seed(12345)
library("nnls")
source("R/signatures.analysis.R")
signatures_decomposition = nmf.fit(contexts_matrix,beta)
save(signatures_decomposition,file="results/signatures_decomposition.RData")
raw_alpha = signatures_decomposition$alpha

# normalize alpha
normalized_alpha = raw_alpha / rowSums(raw_alpha)

# save results
raw_alpha = raw_alpha[,c(3,1,2)]
normalized_alpha = normalized_alpha[,c(3,1,2)]
colnames(raw_alpha) = c("S#1","S#2","S#3")
colnames(normalized_alpha) = c("S#1","S#2","S#3")
save(raw_alpha,file="results/raw_alpha.RData")
save(normalized_alpha,file="results/normalized_alpha.RData")

# make heatmaps
heatmap(normalized_alpha,margin=c(5,5))
