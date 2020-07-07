# load the required libraries
library("factoextra")
library("ggplot2")

# load data and results
load(file="results/normalized_alpha.RData")
load(file="results/contexts_matrix.RData")

# pca analysis to estimate number of clusters
set.seed(11111)
signatures_distance = as.matrix(dist(normalized_alpha))
ydata = prcomp(signatures_distance)
fviz_eig(ydata) + geom_hline(yintercept=5,linetype=3) + geom_text(aes(0,5,label="95% Explained Variance",hjust=-6,vjust=-1)) + geom_hline(yintercept=1,linetype=3) + geom_text(aes(0,1,label="99% Explained Variance",hjust=-6,vjust=-1))

# perform clustering
set.seed(12345)
clusters = t(t(kmeans(ydata$x[,1:3],3,nstart=1000)$cluster))
colnames(clusters) = "Cluster"

# make exposures plot
SUB = clusters[,"Cluster"]
data = NULL
for(i in 1:nrow(normalized_alpha)) {
    curr = cbind(rep(rownames(normalized_alpha)[i],3),as.numeric(normalized_alpha[i,]),names(normalized_alpha[i,]),rep(SUB[i],3))
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

clusters_signatures = t(t(clusters[,"Cluster"]))
rownames(clusters_signatures) = rownames(normalized_alpha)
colnames(clusters_signatures) = "Cluster"
clusters_signatures[which(clusters_signatures[,"Cluster"]==1)] = "Cluster SC#3"
clusters_signatures[which(clusters_signatures[,"Cluster"]==2)] = "Cluster SC#2"
clusters_signatures[which(clusters_signatures[,"Cluster"]==3)] = "Cluster SC#1"
save(clusters_signatures,file="results/clusters_signatures.RData")

curr_data = data.frame(EXPOSURE,SIGNATURE,CLUSTER,stringsAsFactors=FALSE)
ggplot(curr_data,aes(x=SIGNATURE,y=EXPOSURE,fill=CLUSTER)) + geom_boxplot() + xlab("Exposure to Signatures") + ylab("Exposure") + ylim(0,1)
