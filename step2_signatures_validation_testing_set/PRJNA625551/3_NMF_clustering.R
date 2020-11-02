# load the required libraries
library("factoextra")
library("ggplot2")
library("lsa")

# load data and results
load(file="results/raw_alpha.RData")
load(file="results/normalized_alpha.RData")

# perform clustering
set.seed(12345)
clusters = t(t(kmeans(normalized_alpha,3,nstart=1000)$cluster))
colnames(clusters) = "Cluster"

# make exposures plot
SUB = clusters[,"Cluster"]
data = NULL
for(i in 1:nrow(raw_alpha)) {
    curr = cbind(rep(rownames(raw_alpha)[i],3),as.numeric(normalized_alpha[i,]),names(raw_alpha[i,]),rep(SUB[i],3))
    data = rbind(data,curr)
}
colnames(data) = c("PATIENT","EXPOSURE","SIGNATURE","CLUSTER")

EXPOSURE = as.numeric(data[,"EXPOSURE"])
SIGNATURE = as.character(data[,"SIGNATURE"])
SIGNATURE[which(SIGNATURE=="S#1")] = "Signature S#1"
SIGNATURE[which(SIGNATURE=="S#2")] = "Signature S#2"
SIGNATURE[which(SIGNATURE=="S#3")] = "Signature S#3"
CLUSTER = as.character(data[,"CLUSTER"])
CLUSTER[which(CLUSTER==1)] = "Cluster SC#3"
CLUSTER[which(CLUSTER==2)] = "Cluster SC#2"
CLUSTER[which(CLUSTER==3)] = "Cluster SC#1"

signatures_clustering = t(t(clusters[,"Cluster"]))
rownames(signatures_clustering) = rownames(raw_alpha)
colnames(signatures_clustering) = "Cluster"
signatures_clustering[which(signatures_clustering[,"Cluster"]==1)] = "Cluster SC#3"
signatures_clustering[which(signatures_clustering[,"Cluster"]==2)] = "Cluster SC#2"
signatures_clustering[which(signatures_clustering[,"Cluster"]==3)] = "Cluster SC#1"
save(signatures_clustering,file="results/signatures_clustering.RData")

curr_data = data.frame(EXPOSURE,SIGNATURE,CLUSTER,stringsAsFactors=FALSE)
ggplot(curr_data,aes(x=SIGNATURE,y=EXPOSURE,fill=CLUSTER)) + geom_boxplot() + xlab("Exposure to Signatures") + ylab("Exposure") + ylim(0,1)
