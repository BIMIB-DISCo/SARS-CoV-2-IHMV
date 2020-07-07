# load data
load("RData/var_freq_matrix.RData")
load("RData/pvalue_matrix.RData")
load("RData/coverage_matrix.RData")

# set to NA variants with a pvalue >= 0.01
var_freq_matrix[which(pvalue_matrix>=0.01)] = NA

# set to NA variants covered at a coverage <= 10
var_freq_matrix[which(coverage_matrix<=10)] = NA

# select clonal variants
mut_freq = NULL
for(i in 1:ncol(var_freq_matrix)) {
    mut_freq = c(mut_freq,max(var_freq_matrix[,i],na.rm=TRUE))
}
var_freq_matrix = var_freq_matrix[,which(mut_freq>0.90)]

# binarize variants
var_freq_matrix[which(var_freq_matrix<=0.05)] = 0
var_freq_matrix[which(var_freq_matrix>0.90)] = 1
var_freq_matrix[which(var_freq_matrix>0&var_freq_matrix<1)] = NA

# check missing values
missing_values = NULL
for(i in 1:ncol(var_freq_matrix)) {
    missing_values = c(missing_values,length(which(is.na(var_freq_matrix[,i]))))
}
missing_values = missing_values/nrow(var_freq_matrix)

# consider only variants with less than 25% missing values
var_freq_matrix = var_freq_matrix[,which(missing_values<0.25)]

# consider only variants present in at least 5% (approx 50) samples
var_freq_matrix = var_freq_matrix[,names(which((colSums(var_freq_matrix,na.rm=TRUE)/nrow(var_freq_matrix))>0.05))]

# save results
clonal_variants = var_freq_matrix
save(clonal_variants,file="raw_data/clonal_variants.RData")
