# import MrBayes tree
library("ape")
MrBayesTree = read.nexus("results/variants.txt.con.tre")
clades = array(NA,c(length(MrBayesTree$tip.label),1))
rownames(clades) = MrBayesTree$tip.label
colnames(clades) = "Clade"
for(i in 1:MrBayesTree$Nnode) {
    curr_id = i + length(MrBayesTree$tip.label)
    curr_group = MrBayesTree$tip.label[MrBayesTree$edge[which(MrBayesTree$edge[,1]==curr_id),2]]
    clades[curr_group[which(!is.na(curr_group))],"Clade"] = i
}

# compute genotypes for each clade
load("RData/clonal_variants.RData")
samples_thr = 0.03
clonal_variants = clonal_variants[,names(which((colSums(clonal_variants,na.rm=TRUE)/nrow(clonal_variants))>samples_thr))]
genotypes = NULL
for(i in unique(sort(clades[,"Clade"]))) {
    curr_genotypes = clonal_variants[names(which(clades[,"Clade"]==i)),]
    blank_genotype = rep(0,ncol(clonal_variants))
    names(blank_genotype) = colnames(clonal_variants)
    for(j in names(blank_genotype)) {
        blank_genotype[j] = mean(curr_genotypes[,j],na.rm=TRUE)
    }
    genotypes = rbind(genotypes,blank_genotype)
}

# order genotypes based on number of variants
genotypes_ordering = sort.int(rowSums(genotypes),index.return=TRUE,decreasing=FALSE)$ix
genotypes = genotypes[genotypes_ordering,]
rownames(genotypes)[1:9] = paste0("G0",1:9)
rownames(genotypes)[10:25] = paste0("G",10:25)
clades_tmp = clades
cont = 0
for(i in genotypes_ordering) {
    cont = cont + 1
    if(cont<10) {
        curr_group = paste0("G0",cont)
    }
    else {
        curr_group = paste0("G",cont)
    }
    clades[names(which(clades_tmp[,"Clade"]==i)),"Clade"] = curr_group
}

# make oncoprint
library("TRONCO")
data = clonal_variants
data[which(is.na(data))] = 0
data = import.genotypes(data)
data = annotate.stages(data,clades)
oncoprint(data,excl.sort=FALSE,group.by.stage=TRUE)

# save results
save(clades,file="clades.RData")
save(genotypes,file="genotypes.RData")
