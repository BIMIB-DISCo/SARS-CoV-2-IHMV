# load data
library("ape")
MrBayesTree = read.nexus("results/variants.txt.con.tre")
load(file="clades.RData")
load(file="genotypes.RData")

# add clades information to the tree
internal_nodes = rep(NA,MrBayesTree$Nnode)
for(i in sort(unique(clades[,"Clade"]))) {
    curr_tips = which(MrBayesTree$tip.label%in%names(which(clades[,"Clade"]==i)))
    curr_pos = unique(MrBayesTree$edge[which(MrBayesTree$edge[,2]%in%curr_tips),1])-length(MrBayesTree$tip.label)
    internal_nodes[curr_pos] = i
}
MrBayesTree[["node.label"]] = internal_nodes

# re-rooting
edges = MrBayesTree$edge
g01_id = which(MrBayesTree$tip.label%in%names(which(clades[,"Clade"]=="G01")))
g01_id = unique(edges[which(edges[,2]%in%g01_id),1])
MrBayesTree = root(MrBayesTree,node=g01_id)
edges = MrBayesTree$edge
g01_id = which(MrBayesTree$tip.label%in%names(which(clades[,"Clade"]=="G01")))
g01_id = unique(edges[which(edges[,2]%in%g01_id),1]) + 2
innder_ids = (length(MrBayesTree$tip.label)+1):(length(MrBayesTree$tip.label)+MrBayesTree$Nnode)
tips_ids = 1:length(MrBayesTree$tip.label)
edges[which(edges[,1]%in%innder_ids),1] = edges[which(edges[,1]%in%innder_ids),1] + 2
edges[which(edges[,2]%in%innder_ids),2] = edges[which(edges[,2]%in%innder_ids),2] + 2
edges[which(edges[,2]%in%tips_ids),2] = edges[which(edges[,2]%in%tips_ids),2] + 1
root_id = (min(edges[,1])-1)
edge.length = c(c(0.0,0.0),MrBayesTree$edge.length)
edges = rbind(c(root_id,1),edges)
edges = rbind(c(root_id,g01_id),edges)
MrBayesTree$edge = edges
MrBayesTree$edge.length = edge.length
MrBayesTree$Nnode = MrBayesTree$Nnode + 1
MrBayesTree$tip.label = c("SARS-CoV-2-ANC",MrBayesTree$tip.label)
MrBayesTree$node.label = c(MrBayesTree$node.label,"Reference")
MrBayesTree = root(MrBayesTree,node=root_id)

# add clades information to the tree
internal_nodes = rep(NA,MrBayesTree$Nnode)
for(i in sort(unique(clades[,"Clade"]))) {
    curr_tips = which(MrBayesTree$tip.label%in%names(which(clades[,"Clade"]==i)))
    curr_pos = unique(MrBayesTree$edge[which(MrBayesTree$edge[,2]%in%curr_tips),1])-length(MrBayesTree$tip.label)
    internal_nodes[curr_pos] = i
}
MrBayesTree[["node.label"]] = internal_nodes
MrBayesTree$node.label[1] = "Reference"

# write the tree to file in newick format
write.tree(MrBayesTree,file="MrBayesTree.newick",tree.names=TRUE)
