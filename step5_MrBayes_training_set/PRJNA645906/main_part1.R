# settings
variants_thr_low = 0.05
variants_thr_high = 0.90
samples_thr = 0.03

# load variants data
load("RData/clonal_variants.RData")
load("RData/processed_variants.RData")

# select only frequent clonal variants
clonal_variants = clonal_variants[,names(which((colSums(clonal_variants,na.rm=TRUE)/nrow(clonal_variants))>samples_thr))]

# build MrBayes input file
variants = processed_variants[,colnames(clonal_variants)]
variants = variants[rownames(clonal_variants),]
present = which(variants>variants_thr_high,arr.ind=TRUE)
absent = which(variants<variants_thr_low,arr.ind=TRUE)
ambiguous = which(variants>=variants_thr_low&variants<=variants_thr_high,arr.ind=TRUE)
variants[present] = "1"
variants[absent] = "0"
variants[ambiguous] = "{0,1}"
variants[which(is.na(variants))] = "{0,1}"

# write results to file
write.table(x=variants,file="variants.txt",quote=FALSE,sep="\t",row.names=TRUE,col.names=FALSE)
