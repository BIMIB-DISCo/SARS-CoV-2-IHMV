# load required libraries and sources
library("ape")
library("igraph")
source("R/VERSOphylo.R")

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

# compute phylogeny by VERSO
adj_matrix = as.adj.matrix(inference$B[-1,-1])
phylogeny = VERSOphylo(adj_matrix,corrected_genotype[,rownames(adj_matrix)])

# write the tree to file in newick format
write.tree(phylogeny,file="VERSO.newick",tree.names=TRUE)
