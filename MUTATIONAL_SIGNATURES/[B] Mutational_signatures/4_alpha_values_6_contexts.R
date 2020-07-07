# load results
load("results/contexts_matrix.RData")
load("results/signatures_results.RData")
observed_counts = contexts_matrix
raw_alpha = signatures_results$decomposition$alpha[[3]]
predicted_counts = signatures_results$decomposition$alpha[[3]] %*% signatures_results$decomposition$beta[[3]]

# consider only samples at high correlation between observed and predicted counts (correlation >0.90)
raw_alpha = raw_alpha[names(which(diag(cor(t(observed_counts),t(predicted_counts)))>0.90)),]
print(round(nrow(raw_alpha)/nrow(observed_counts),digits=2)) # keep approx 96% of patients

# normalize alpha
normalized_alpha = raw_alpha / rowSums(raw_alpha)

# save results
raw_alpha = raw_alpha[,c(2,1,3)]
normalized_alpha = normalized_alpha[,c(2,1,3)]
colnames(raw_alpha) = c("S#1","S#2","S#3")
colnames(normalized_alpha) = c("S#1","S#2","S#3")
save(raw_alpha,file="results/raw_alpha.RData")
save(normalized_alpha,file="results/normalized_alpha.RData")

# make heatmaps
heatmap(normalized_alpha,margin=c(5,5))
