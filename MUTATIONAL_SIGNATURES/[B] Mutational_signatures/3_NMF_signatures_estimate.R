# load the required libraries
library("ggplot2")
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

# load data and results
load(file="raw_data/reference.RData")
load(file="results/contexts_matrix.RData")
load(file="results/signatures_results.RData")
load(file="results/trinucleotides_matrix.RData")

# number of signatures

RANGE = 1:6
VALUE = as.numeric(signatures_results$decomposition$measures[,"Cophenetic_Coefficient"])
VALUE[1] = 1
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("Evaluation of Optimal Number of Signatures") + xlab("Number of Signatures") + ylab("Cophenetic Coefficient (1000 NMF restarts)")

RANGE = 1:6
VALUE = as.numeric(signatures_results$decomposition$measures[,"Dispersion_Coefficient"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("Evaluation of Optimal Number of Signatures") + xlab("Number of Signatures") + ylab("Dispersion Coefficient (1000 NMF restarts)")

RANGE = 1:6
VALUE = as.numeric(signatures_results$decomposition$measures[,"Explained_Variance"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("Evaluation of Optimal Number of Signatures") + xlab("Number of Signatures") + ylab("Explained Variance (1000 NMF restarts)")

RANGE = 1:6
VALUE = NULL
observations = contexts_matrix
for(i in 1:6) {
    predictions = signatures_results$decomposition$alpha[[i]]%*%signatures_results$decomposition$beta[[i]]
    VALUE = c(VALUE,mean(diag(cor(t(observations),t(predictions)))))
}
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("Evaluation of Optimal Number of Signatures") + xlab("Number of Signatures") + ylab("Goodness of fit")

RANGE = 1:6
VALUE = as.numeric(signatures_results$cross_validation$summary[,"Mean cross-validation error"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("Evaluation of Optimal Number of Signatures") + xlab("Number of Signatures") + ylab("Mean cross-validation error (1000 cross-validation Iterations)")

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

# 6 contexts signatures
beta = signatures_results$decomposition$beta[[3]]
colnames(beta) = c("C>A | G>T","C>G | G>C","C>T | G>A","T>A | A>T","T>C | A>G","T>G | A>C")
data = rbind(theoretical_6_contexts,beta)
data = data[c(1,3,2,4),]
rownames(data) = c("Reference","Signature S#1","Signature S#2","Signature S#3")
plot.6.contexts(data)
