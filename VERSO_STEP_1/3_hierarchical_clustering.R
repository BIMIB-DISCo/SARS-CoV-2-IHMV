# load data
load(file="RData/samples_similarity.RData")

# get clusters
row.clusters = hclust(as.dist((max(samples_similarity)-samples_similarity)))
assignments = cutree(row.clusters,k=11)
assignments = t(t(assignments))
colnames(assignments) = "Cluster"

# relabeling
assignments[which(assignments[,1]==4),1] = "C01"
assignments[which(assignments[,1]==3),1] = "C02"
assignments[which(assignments[,1]==1),1] = "C03"
assignments[which(assignments[,1]==2),1] = "C04"
assignments[which(assignments[,1]==9),1] = "C05"
assignments[which(assignments[,1]==8),1] = "C06"
assignments[which(assignments[,1]==6),1] = "C07"
assignments[which(assignments[,1]==11),1] = "C08"
assignments[which(assignments[,1]==5),1] = "C09"
assignments[which(assignments[,1]==10),1] = "C10"
assignments[which(assignments[,1]==7),1] = "C11"
save(assignments,file="assignments.RData")
