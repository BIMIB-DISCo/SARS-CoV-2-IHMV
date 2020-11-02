# load data
load("results/var_freq_matrix.RData")
load("results/alt_count_matrix.RData")
load("results/pvalue_matrix.RData")
load("results/coverage_matrix.RData")

# set to 0 variants with frequency <= 0.05
var_freq_matrix[which(var_freq_matrix<=0.05)] = 0

# set to 0 variants with alt count < 3
var_freq_matrix[which(alt_count_matrix<3)] = 0

# set to NA variants with pvalue >= 0.01
var_freq_matrix[which(pvalue_matrix>=0.01)] = NA

# set to NA variants with coverage < 20
var_freq_matrix[which(coverage_matrix<20)] = NA

# check missing values
missing_values = NULL
for(i in 1:ncol(var_freq_matrix)) {
    missing_values = c(missing_values,length(which(is.na(var_freq_matrix[,i]))))
}
missing_values = missing_values/nrow(var_freq_matrix)

# consider only variants with less than 25% missing values
var_freq_matrix = var_freq_matrix[,which(missing_values<0.25)]

# consider only variants present in at least 1 sample
mutated_variants = NULL
for(i in 1:ncol(var_freq_matrix)) {
    mutated_variants = c(mutated_variants,length(which(var_freq_matrix[,i]>0)))
}
var_freq_matrix = var_freq_matrix[,which((mutated_variants/nrow(var_freq_matrix))>0.00)]

# save results
processed_variants = var_freq_matrix
save(processed_variants,file="processed_variants.RData")
