# load the libraries
library("gplots")
library("RColorBrewer")
library("cluster")
library("factoextra")
library("ggplot2")
library("igraph")
library("SparseSignatures")
library("data.table")
library("gridExtra")
"plot.6.contexts" = function(beta) {
    x <- as.data.table(melt(beta, varnames = c("signature", "cat")))
    glist <- list()
    for (i in 1:nrow(beta)) {
        plt <- ggplot(x[signature == rownames(beta)[i]]) + geom_bar(aes(x = cat, 
            y = value, fill = cat), stat = "identity", position = "identity") + 
            facet_wrap(~cat, nrow = 1, scales = "free_x") + theme(axis.text.x = element_text(angle = 90, 
            hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
            ggtitle(rownames(beta)[i]) + theme(text = element_text(size=15),legend.position = "none") + 
            ylab("Frequency of mutations") + xlab("") + ylim(0,1)
        plt <- plt + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        glist[[i]] <- plt
    }
    grid.arrange(grobs = glist, ncol = ceiling(nrow(beta)/3))
}
"plot.96.contexts" = function( beta, useColNames = TRUE, mutation_categories = NULL) 
{
    x <- as.data.table(melt(beta, varnames = c("signature", "cat")))
    x[, `:=`(Context, paste0(substr(cat, 1, 1), ".", substr(cat, 
        7, 7)))]
    x[, `:=`(alt, paste0(substr(cat, 3, 3), ">", substr(cat, 
        5, 5)))]
    x$alt[which(x$alt=="C>A")] = "C>A | G>T"
    x$alt[which(x$alt=="C>G")] = "C>G | G>C"
    x$alt[which(x$alt=="C>T")] = "C>T | G>A"
    x$alt[which(x$alt=="T>A")] = "T>A | A>T"
    x$alt[which(x$alt=="T>C")] = "T>C | A>G"
    x$alt[which(x$alt=="T>G")] = "T>G | A>C"
    glist <- list()
    for (i in 1:nrow(beta)) {
        plt <- ggplot(x[signature == rownames(beta)[i]]) + geom_bar(aes(x = Context, 
            y = value, fill = alt), stat = "identity", position = "identity") + 
            facet_wrap(~alt, nrow = 1, scales = "free_x") + theme(axis.text.x = element_text(angle = 90, 
            hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
            ggtitle(rownames(beta)[i]) + theme(text = element_text(size=15),legend.position = "none") + 
            ylab("Frequency of mutations") + ylim(0,0.10)
        plt <- plt + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
        glist[[i]] <- plt
    }
    grid.arrange(grobs = glist, ncol = ceiling(nrow(beta)/3))
}

# load data and results
load(file="processed_data/reference.RData")
load(file="results/contexts_matrix.RData")
load(file="results/trinucleotides_matrix.RData")
load(file="results/normalized_alpha.RData")
load(file="results/signatures_clustering.RData")

# make heatmap for contexts
my_palette = colorRampPalette(c("blue","white","red"))(n=300)
mat_data = normalized_alpha[rownames(signatures_clustering),]
mat_data = mat_data[names(sort(gsub("Cluster SC#","",signatures_clustering[,"Cluster"]))),]
clustering = as.numeric(sort(gsub("Cluster SC#","",signatures_clustering[,"Cluster"])))
tmp = clustering
clustering[which(tmp==1)] = 1
clustering[which(tmp==2)] = 2
clustering[which(tmp==3)] = 3
mat_data = data.matrix(mat_data[sort.int(clustering,index.return=TRUE)$ix,])
clustering = sort(clustering)
heatmap.2(t(mat_data),main="Mutational Signatures Clusters",notecol="black",density.info="none",trace="none",margins=c(12,9),col=my_palette,dendrogram="row",Rowv=TRUE,Colv=FALSE,ColSideColors=as.character(clustering))
legend("topright",legend=c("Cluster SC#1","Cluster SC#2","Cluster SC#3"),col=unique(clustering),lty=1,lwd=5,cex=.7)

# compute theoretical trinucleotides
trinucleotides = colnames(trinucleotides_matrix)
theoretical_trinucleotides = array(0,c(1,length(trinucleotides)))
rownames(theoretical_trinucleotides) = "Reference"
colnames(theoretical_trinucleotides) = trinucleotides
for(i in 2:(length(reference)-1)) {
    valid_muts = c("A","C","G","T")[which(!c("A","C","G","T")==reference[i])]
    for(j in valid_muts) {
        curr_new_mut = paste0(reference[(i-1)],"[",reference[i],">",j,"]",reference[(i+1)])
        # merge equivalent contexts
        curr_new_mut = gsub("G>T","C>A",curr_new_mut)
        curr_new_mut = gsub("G>C","C>G",curr_new_mut)
        curr_new_mut = gsub("G>A","C>T",curr_new_mut)
        curr_new_mut = gsub("A>T","T>A",curr_new_mut)
        curr_new_mut = gsub("A>G","T>C",curr_new_mut)
        curr_new_mut = gsub("A>C","T>G",curr_new_mut)
        theoretical_trinucleotides["Reference",curr_new_mut] = theoretical_trinucleotides["Reference",curr_new_mut] + 1
    }
}
theoretical_96_contexts = theoretical_trinucleotides / rowSums(theoretical_trinucleotides)
theoretical_6_contexts = array(NA,c(1,6))
rownames(theoretical_6_contexts) = "Reference"
colnames(theoretical_6_contexts) = c("C>A | G>T","C>G | G>C","C>T | G>A","T>A | A>T","T>C | A>G","T>G | A>C")
theoretical_6_contexts["Reference","C>A | G>T"] = sum(theoretical_trinucleotides[,grep("C>A",colnames(theoretical_trinucleotides))])
theoretical_6_contexts["Reference","C>G | G>C"] = sum(theoretical_trinucleotides[,grep("C>G",colnames(theoretical_trinucleotides))])
theoretical_6_contexts["Reference","C>T | G>A"] = sum(theoretical_trinucleotides[,grep("C>T",colnames(theoretical_trinucleotides))])
theoretical_6_contexts["Reference","T>A | A>T"] = sum(theoretical_trinucleotides[,grep("T>A",colnames(theoretical_trinucleotides))])
theoretical_6_contexts["Reference","T>C | A>G"] = sum(theoretical_trinucleotides[,grep("T>C",colnames(theoretical_trinucleotides))])
theoretical_6_contexts["Reference","T>G | A>C"] = sum(theoretical_trinucleotides[,grep("T>G",colnames(theoretical_trinucleotides))])
theoretical_6_contexts = theoretical_6_contexts / rowSums(theoretical_6_contexts)

# 6 contexts clusters
data = rbind(colSums(contexts_matrix[names(which(gsub("Cluster SC#","",signatures_clustering[,"Cluster"])=="1")),]),colSums(contexts_matrix[names(which(gsub("Cluster SC#","",signatures_clustering[,"Cluster"])=="2")),]),colSums(contexts_matrix[names(which(gsub("Cluster SC#","",signatures_clustering[,"Cluster"])=="3")),]))
data = rbind(theoretical_6_contexts,data)
data = data/rowSums(data)
rownames(data) = c("Reference","Cluster SC#1","Cluster SC#2","Cluster SC#3")
plot.6.contexts(data)

# 96 trinucleotides clusters
data = rbind(colSums(trinucleotides_matrix[names(which(gsub("Cluster SC#","",signatures_clustering[,"Cluster"])=="1")),]),colSums(trinucleotides_matrix[names(which(gsub("Cluster SC#","",signatures_clustering[,"Cluster"])=="2")),]),colSums(trinucleotides_matrix[names(which(gsub("Cluster SC#","",signatures_clustering[,"Cluster"])=="3")),]))
data = rbind(theoretical_96_contexts,data)
data = data/rowSums(data)
rownames(data) = c("Reference","Cluster SC#1","Cluster SC#2","Cluster SC#3")
plot.96.contexts(data)
