# load results
load("results/contexts_matrix.RData")
load("results/signatures_decomposition.RData")
observed_counts = contexts_matrix
raw_alpha = signatures_decomposition$alpha[[3]]

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
