# load the required libraries
library("factoextra")
library("TRONCO")

# load data
load(file="raw_data/clonal_variants.RData")
load(file="results/best_LACE.RData")
load(file="results/C_extended.RData")
load(file="results/samples_similarity.RData")

# save observed and corrected genotypes
observed_genotype = clonal_variants[rownames(C_extended),]
B = best_LACE$B[-1,-1]
C = C_extended[,colnames(B)]
corrected_genotype = C%*%B
corrected_genotype = corrected_genotype[rownames(observed_genotype),]
corrected_genotype = corrected_genotype[,colnames(observed_genotype)]

# compute errors estimate (FP and FN)
errors_estimate = array(NA,c(ncol(observed_genotype),5))
rownames(errors_estimate) = colnames(observed_genotype)
colnames(errors_estimate) = c("True_Positive","True_Negative","False_Positive","False_Negative","Hamming_Distance")
for(i in 1:nrow(errors_estimate)) {
    curr_observed = as.numeric(observed_genotype[,rownames(errors_estimate)[i]])
    curr_corrected = as.numeric(corrected_genotype[,rownames(errors_estimate)[i]])
    errors_estimate[i,"True_Positive"] = length(which(curr_observed==1&curr_corrected==1))
    errors_estimate[i,"True_Negative"] = length(which(curr_observed==0&curr_corrected==0))
    errors_estimate[i,"False_Positive"] = length(which(curr_observed==1&curr_corrected==0))
    errors_estimate[i,"False_Negative"] = length(which(curr_observed==0&curr_corrected==1))
    errors_estimate[i,"Hamming_Distance"] = errors_estimate[i,"False_Positive"] + errors_estimate[i,"False_Negative"]
}
print(sort(errors_estimate[,"False_Positive"])/ncol(observed_genotype))
print(sort(errors_estimate[,"False_Negative"])/ncol(observed_genotype))
print(sort(errors_estimate[,"Hamming_Distance"])/ncol(observed_genotype))

# pca analysis
set.seed(12345)
ydata = prcomp(as.dist((max(samples_similarity)-samples_similarity)))
fviz_eig(ydata,ncp=24)

# get clusters
set.seed(54321)
row.clusters = hclust(as.dist((max(samples_similarity)-samples_similarity)))
assignments = cutree(row.clusters,k=13)
assignments = t(t(assignments))
colnames(assignments) = "Cluster"

# relabeling
assignments[which(assignments[,1]==1),1] = "C13"
assignments[which(assignments[,1]==2),1] = "C09"
assignments[which(assignments[,1]==3),1] = "C06"
assignments[which(assignments[,1]==4),1] = "C01"
assignments[which(assignments[,1]==5),1] = "C03"
assignments[which(assignments[,1]==6),1] = "C12"
assignments[which(assignments[,1]==7),1] = "C10"
assignments[which(assignments[,1]==8),1] = "C04"
assignments[which(assignments[,1]==9),1] = "C05"
assignments[which(assignments[,1]==10),1] = "C07"
assignments[which(assignments[,1]==11),1] = "C02"
assignments[which(assignments[,1]==12),1] = "C11"
assignments[which(assignments[,1]==13),1] = "C08"
save(assignments,file="results/assignments.RData")

# oncoprint NA to 0
load(file="raw_data/clonal_variants.RData")
clonal_variants = clonal_variants[rownames(assignments),]
clonal_variants[which(is.na(clonal_variants))] = 0
data = import.genotypes(clonal_variants)
data = annotate.stages(data,assignments)
oncoprint(data,excl.sort=FALSE,group.by.stage=TRUE)

# oncoprint NA to 1
load(file="raw_data/clonal_variants.RData")
clonal_variants = clonal_variants[rownames(assignments),]
clonal_variants[which(is.na(clonal_variants))] = 1
data = import.genotypes(clonal_variants)
data = annotate.stages(data,assignments)
oncoprint(data,excl.sort=FALSE,group.by.stage=TRUE)
