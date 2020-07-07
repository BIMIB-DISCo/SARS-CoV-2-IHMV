# load required libraries and scripts
library("parallel")
library("Rfast")
library("igraph")
source("R/frontend.R")
source("R/inference.R")
source("R/utils.R")
source("R/visualization.R")
"get.adj.matrix" <- function( B ) {
    
    adj_matrix <- array(0L,c((dim(B)[1]-1),(dim(B)[2]-1)))
    rownames(adj_matrix) <- colnames(B)[2:ncol(B)]
    colnames(adj_matrix) <- colnames(B)[2:ncol(B)]
    
    for(rP in 2:(nrow(B)-1)) {
        
        for(rC in ((rP+1):nrow(B))) {
            
            if(all(B[rP,1:rP]==B[rC,1:rP])&&(sum(B[rP,])==(sum(B[rC,])-1))) {
                
                adj_matrix[(rP-1),(rC-1)] <- 1
                
            }
            
        }
        
    }
    
    return(adj_matrix)
    
}

# load results
load(file="results/inference.RData")

# consensus matrix
cont = 0
B = inference$B
consensus_matrix = get.adj.matrix(B)
cont = cont + 1
for(i in 1:length(inference$equivalent_solutions)) {
    curr_adj_matrix = get.adj.matrix(inference$equivalent_solutions[[i]][["B"]])
    curr_adj_matrix = curr_adj_matrix[rownames(consensus_matrix),colnames(consensus_matrix)]
    consensus_matrix = consensus_matrix + curr_adj_matrix
    cont = cont + 1
}
consensus_matrix = consensus_matrix / cont

# score all maximum likelihood trees
ml_trees = NULL
cont = 0
B = inference$B
curr_adj_matrix = get.adj.matrix(B)
curr_adj_matrix = curr_adj_matrix[rownames(consensus_matrix),colnames(consensus_matrix)]
curr_adj_matrix = curr_adj_matrix * consensus_matrix
cont = cont + 1
ml_trees = c(ml_trees,(sum(curr_adj_matrix)/(nrow(curr_adj_matrix)-1)))
for(i in 1:length(inference$equivalent_solutions)) {
    curr_adj_matrix = get.adj.matrix(inference$equivalent_solutions[[i]][["B"]])
    curr_adj_matrix = curr_adj_matrix[rownames(consensus_matrix),colnames(consensus_matrix)]
    curr_adj_matrix = curr_adj_matrix * consensus_matrix
    cont = cont + 1
    ml_trees = c(ml_trees,(sum(curr_adj_matrix)/(nrow(curr_adj_matrix)-1)))
}

# make plots

dev.new(height=15,width=15)
heatmap(consensus_matrix)

dev.new(height=15,width=15)
best_LACE = inference$equivalent_solutions[[(which.max(ml_trees)-1)]]
inference$B = best_LACE$B
inference$C = best_LACE$C
best_LACE = inference
longitudinal.tree.plot(best_LACE,show_prev=FALSE)

# save results
save(consensus_matrix,file="results/consensus_matrix.RData")
save(best_LACE,file="results/best_LACE.RData")
