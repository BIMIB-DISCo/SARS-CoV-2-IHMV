# load required libraries and sources
library("ape")
library("igraph")
library("TRONCO")
source("R/VERSOphylo.R")

# set the seed
set.seed(33333)

# load data and results
load(file="results/clonal_variants.RData")
clonal_variants = clonal_variants[,-which(colnames(clonal_variants)=="28883_G_C")]
colnames(clonal_variants)[which(colnames(clonal_variants)=="28882_G_A")] = "28882_G_A|28883_G_C"
load(file="results/inference_step2.RData")
inference = inference_step2

# compute observed and corrected genotypes
observed_genotype = clonal_variants[,colnames(inference$B[-1,-1])]
B = inference$B[-1,-1]
C_vector = inference$C
C = array(0,c(nrow(clonal_variants),ncol(B)))
rownames(C) = rownames(clonal_variants)
colnames(C) = colnames(B)
for(i in 1:nrow(C)) {
    C[i,(as.numeric(gsub("G","",C_vector[i]))-1)] = 1
}
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

# compute phylogeny by VERSO
adj_matrix = as.adj.matrix(inference$B[-1,-1])
phylogeny = VERSOphylo(adj_matrix,corrected_genotype[,rownames(adj_matrix)])

# compute cluster assignments
phylogeny = drop.tip(phylogeny,phylogeny$tip.label[which(phylogeny$tip.label=="Reference"):length(phylogeny$tip.label)])
distance = cophenetic.phylo(phylogeny)
row.clusters = hclust(as.dist(distance),method="complete")
assignments = cutree(row.clusters,k=length(unique(inference$C[,1])))
assignments = t(t(assignments))
colnames(assignments) = "Cluster"

# determine ordering of clusters
num_variants = NULL
for(i in sort(unique(assignments[,"Cluster"]))) {
    num_variants = c(num_variants,sum(colSums(clonal_variants[names(which(assignments[,"Cluster"]==i)),,drop=FALSE],na.rm=TRUE)/length(names(which(assignments[,"Cluster"]==i)))))
}
genotypes_ordering = sort.int(num_variants,index.return=TRUE)$ix

# relabeling
cont = 0
for(i in genotypes_ordering) {
    cont = cont + 1
    if(cont<10) {
        curr_cluster = paste0("C0",cont)
    }
    else {
        curr_cluster = paste0("C",cont)
    }
    assignments[which(assignments[,1]==i),1] = curr_cluster
}
save(assignments,file="results/assignments.RData")

# oncoprint NA to 0
data = clonal_variants[,colnames(observed_genotype)]
data = data[rownames(observed_genotype),]
data[which(is.na(data))] = 0
data = import.genotypes(data)
data = annotate.stages(data,assignments)
oncoprint(data,excl.sort=FALSE,group.by.stage=TRUE)

# oncoprint NA to 1
data = clonal_variants[,colnames(observed_genotype)]
data = data[rownames(observed_genotype),]
data[which(is.na(data))] = 1
data = import.genotypes(data)
data = annotate.stages(data,assignments)
oncoprint(data,excl.sort=FALSE,group.by.stage=TRUE)
