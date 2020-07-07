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
"compute.likelihood.attachament" <- function( D, B, curr_alpha, curr_beta ) {

    # Go through all the cells
    D[which(is.na(D))] = -3
    curr_cells_D <- D
    
    # Find assignment at maximum log likelihood for current cell
    lik_matrix <- array(0L,c(nrow(D),(ncol(B)-1)))
    
    for(k in 2:nrow(B)) {
        
        curr_clone_C = matrix(rep(0L,(nrow(B)-1)),nrow=1)
        curr_clone_C[1,(k-1)] <- 1L
        
        # Save mapping between ordering in B and dataset D
        # Factorization D = C dot B. r_D_tilde represents a combination of mutations and clones
        r_D_tilde <- (curr_clone_C %*% B[-1,-1,drop=FALSE])*2
        
        sum_cell_clone <- as.matrix(sweep(curr_cells_D,MARGIN=2,r_D_tilde,"+"))
        lik_matrix[,(k-1)] <- (curr_beta^rowSums(sum_cell_clone==2)) * ((1-curr_beta)^rowSums(sum_cell_clone==0)) * ((curr_alpha)^rowSums(sum_cell_clone==1)) * ((1-curr_alpha)^rowSums(sum_cell_clone==3))
        
    }

    return(lik_matrix)

}

# load data and results
load(file="raw_data/clonal_variants.RData")
load(file="results/inference.RData")

# consensus matrix
samples_similarity = array(0,c(nrow(inference$C$Experiment_1),nrow(inference$C$Experiment_1)))
rownames(samples_similarity) = rownames(inference$C$Experiment_1)
colnames(samples_similarity) = rownames(inference$C$Experiment_1)

# consider all inferred C
cont = 0
B = inference$B
C = inference$C[[1]]
C_extended = array(0,c(nrow(C),(nrow(B)-1)))
rownames(C_extended) = rownames(C)
colnames(C_extended) = colnames(B)[-1]
likelihood = compute.likelihood.attachament(clonal_variants[,colnames(B)[-1]],B,inference$error_rates$alpha,inference$error_rates$beta)
for(i in 1:nrow(C_extended)) {
    C_extended[i,colnames(B)[-1][which(likelihood[i,]==max(likelihood[i,]))]] = C_extended[i,colnames(B)[-1][which(likelihood[i,]==max(likelihood[i,]))]] + 1
}
for(a in rownames(samples_similarity)) {
    for(b in colnames(samples_similarity)) {
        samples_similarity[a,b] = samples_similarity[a,b] + length(which(C_extended[a,]==1&C_extended[b,]==1))
    }
}
C_extended_full = C_extended
cont = cont + 1
for(i in 1:length(inference$equivalent_solutions)) {
    B = inference$equivalent_solutions[[i]]$B
    C = matrix(inference$equivalent_solutions[[i]]$C[[1]],ncol=1)
    rownames(C) = rownames(C_extended_full)
    colnames(C) = "Clone"
    C_extended = array(0,c(nrow(C),(nrow(B)-1)))
    rownames(C_extended) = rownames(C)
    colnames(C_extended) = colnames(B)[-1]
    likelihood = compute.likelihood.attachament(clonal_variants[,colnames(B)[-1]],B,inference$error_rates$alpha,inference$error_rates$beta)
    for(i in 1:nrow(C_extended)) {
        C_extended[i,colnames(B)[-1][which(likelihood[i,]==max(likelihood[i,]))]] = C_extended[i,colnames(B)[-1][which(likelihood[i,]==max(likelihood[i,]))]] + 1
    }
    for(a in rownames(samples_similarity)) {
        for(b in colnames(samples_similarity)) {
            samples_similarity[a,b] = samples_similarity[a,b] + length(which(C_extended[a,]==1&C_extended[b,]==1))
        }
    }
    C_extended_full = C_extended_full + C_extended[,colnames(C_extended_full)]
    cont = cont + 1
}
C_extended = C_extended_full

# make plots
dev.new(height=15,width=15)
heatmap(C_extended)
dev.new(height=15,width=15)
heatmap(samples_similarity)

# save results
save(C_extended,file="results/C_extended.RData")
save(samples_similarity,file="results/samples_similarity.RData")
