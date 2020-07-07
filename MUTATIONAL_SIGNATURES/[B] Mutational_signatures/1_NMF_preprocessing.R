# load data
load(file="raw_data/processed_variants.RData")
load(file="raw_data/reference.RData")
load(file="raw_data/trinucleotides.RData")

# remove clonal variants with VF>0.90 in at least one sample
mut_freq = NULL
for(i in 1:ncol(processed_variants)) {
    mut_freq = c(mut_freq,max(processed_variants[,i],na.rm=TRUE))
}
valid_variants = processed_variants[,which(mut_freq<=0.90)]
valid_variants[which(valid_variants<=0.05)] = 0

# select patients and variants contributing at least 6 mutations
binarized_valid_variants = valid_variants
binarized_valid_variants[which(is.na(binarized_valid_variants))] = 0
binarized_valid_variants[which(binarized_valid_variants>0)] = 1
valid_variants = valid_variants[names(which(rowSums(binarized_valid_variants)>=6)),]
valid_variants = valid_variants[,names(which(colSums(valid_variants,na.rm=TRUE)>0.00))]

# get trinucleotides contexts
trinucleotide_variants = NULL
for(v in colnames(valid_variants)) {
    # get trinucleotide context of current variant
    curr_split = strsplit(v,split="_")[[1]]
    curr_trinucleotide = paste0(reference[(as.numeric(curr_split[1])-1)],"[",curr_split[2],">",curr_split[3],"]",reference[(as.numeric(curr_split[1])+1)])
    trinucleotide_variants = c(trinucleotide_variants,curr_trinucleotide)
}
colnames(valid_variants) = trinucleotide_variants

# built 96 trinucleotides matrix
trinucleotides_matrix = array(NA,c(nrow(valid_variants),96))
rownames(trinucleotides_matrix) = rownames(valid_variants)
colnames(trinucleotides_matrix) = trinucleotides
valid_variants_trinucleotides = valid_variants
colnames(valid_variants_trinucleotides) = gsub("G>T","C>A",colnames(valid_variants_trinucleotides))
colnames(valid_variants_trinucleotides) = gsub("G>C","C>G",colnames(valid_variants_trinucleotides))
colnames(valid_variants_trinucleotides) = gsub("G>A","C>T",colnames(valid_variants_trinucleotides))
colnames(valid_variants_trinucleotides) = gsub("A>T","T>A",colnames(valid_variants_trinucleotides))
colnames(valid_variants_trinucleotides) = gsub("A>G","T>C",colnames(valid_variants_trinucleotides))
colnames(valid_variants_trinucleotides) = gsub("A>C","T>G",colnames(valid_variants_trinucleotides))
for(i in 1:nrow(trinucleotides_matrix)) {
    for(j in colnames(trinucleotides_matrix)) {
        trinucleotides_matrix[i,j] = sum(valid_variants_trinucleotides[i,which(colnames(valid_variants_trinucleotides)==j)],na.rm=TRUE)
    }
}

# built 6 contexts matrix
contexts_matrix = array(NA,c(nrow(trinucleotides_matrix),6))
rownames(contexts_matrix) = rownames(trinucleotides_matrix)
colnames(contexts_matrix) = c("C>A","C>G","C>T","T>A","T>C","T>G")
for(i in 1:nrow(contexts_matrix)) {
    for(j in colnames(contexts_matrix)) {
        contexts_matrix[i,j] = sum(trinucleotides_matrix[i,grep(j,colnames(trinucleotides_matrix))],na.rm=TRUE)
    }
}

# select patients with a cumulative VF greater than 1
valid_patients = sort(unique(names(which(rowSums(contexts_matrix)>1))))
contexts_matrix = contexts_matrix[valid_patients,]
trinucleotides_matrix = trinucleotides_matrix[valid_patients,]

# save results
save(contexts_matrix,file="results/contexts_matrix.RData")
save(trinucleotides_matrix,file="results/trinucleotides_matrix.RData")
