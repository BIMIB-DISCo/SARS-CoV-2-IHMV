VERSO_v2 <- function( D, alpha = NULL, beta = NULL, check_indistinguishable = TRUE, num_rs = 10, num_iter = 10000, n_try_bs = 1000, num_processes = Inf, verbose = TRUE, initialization = NULL ) {
    
    # set storage mode to integer
    storage.mode(D) <- "integer"
    
    # remove any indistinguishable variant from input data prior inference
    if(check_indistinguishable) {
        D <- check.indistinguishable(D)
        storage.mode(D) <- "integer"
    }

    # initialize error rates alpha and beta if not provided as inputs
    if(is.null(alpha)) {
        alpha = c(0.01,0.05,0.10)
        beta = c(0.01,0.05,0.10)
    }
    
    if(verbose) {
        message("Performing inference for a total of ",length(alpha)," different values of alpha and beta.")
    }

    # setting up parallel execution
    parallel <- NULL
    close_parallel <- FALSE
    if(is.null(parallel)&&length(alpha)>1) {
        if(is.na(num_processes)||is.null(num_processes)||num_processes==1) {
            parallel <- NULL
        }
        else if(num_processes==Inf) {
            cores <- as.integer((detectCores()-1))
            if(cores < 2) {
                parallel <- NULL
            }
            else {
                num_processes <- min(num_processes,length(alpha))
                parallel <- makeCluster(num_processes,outfile="")
                close_parallel <- TRUE
            }
        }
        else {
            num_processes <- min(num_processes,length(alpha))
            parallel <- makeCluster(num_processes,outfile="")
            close_parallel <- TRUE
        }
        if(verbose && !is.null(parallel)) {
            message("Executing ",num_processes," processes in parallel.")
        }
    }

    # set initial tree from where to start MCMC search
    # data <- D
    # data[which(is.na(data))] <- 0
    # marginal_probs <- matrix(colSums(data,na.rm=TRUE)/nrow(data),ncol=1)
    # rownames(marginal_probs) <- colnames(data)
    # colnames(marginal_probs) <- "Frequency"
    # joint_probs <- array(NA,c(ncol(data),ncol(data)))
    # rownames(joint_probs) <- colnames(data)
    # colnames(joint_probs) <- colnames(data)
    # for (i in seq_len(ncol(data))) {
    #     for (j in seq_len(ncol(data))) {
    #         val1 <- data[,i]
    #         val2 <- data[,j]
    #         joint_probs[i,j] <- (t(val1)%*%val2)
    #     }
    # }
    # joint_probs <- joint_probs/nrow(data)
    # adjacency_matrix <- array(0,c(ncol(data),ncol(data)))
    # rownames(adjacency_matrix) <- colnames(data)
    # colnames(adjacency_matrix) <- colnames(data)
    # pmi <- joint_probs
    # for(i in seq_len(nrow(pmi))) {
    #     for(j in seq_len(ncol(pmi))) {
    #         pmi[i,j] <- log(joint_probs[i,j]/(marginal_probs[i,"Frequency"]*marginal_probs[j,"Frequency"]))
    #     }
    # }
    # ordering <- names(sort(marginal_probs[,"Frequency"],decreasing=TRUE))
    # adjacency_matrix <- adjacency_matrix[ordering,ordering]
    # adjacency_matrix[1,2] = 1
    # if(nrow(adjacency_matrix)>2) {
    #     for(i in 3:nrow(adjacency_matrix)) {
    #         curr_c <- rownames(adjacency_matrix)[i]
    #         curr_candidate_p <- rownames(adjacency_matrix)[seq_len((i-1))]
    #         adjacency_matrix[names(which.max(pmi[curr_candidate_p,curr_c]))[1],curr_c] <- 1
    #     }
    # }
    # adjacency_matrix <- rbind(rep(0,nrow(adjacency_matrix)),adjacency_matrix)
    # adjacency_matrix <- cbind(rep(0,nrow(adjacency_matrix)),adjacency_matrix)
    # adjacency_matrix[1,2] = 1
    # initialization <- as.B(adj_matrix=adjacency_matrix,D=D)

    # now start the inference
    if(is.null(parallel)) {

        # sequential computation
        inference <- list()
        for(i in seq_len(length(alpha))) {
            
            if(verbose) {
                message(paste('Performing inference for alpha =',paste0(alpha[i],collapse=" | "),'and beta =',paste0(beta[i],collapse=" | ")))
            }

            inference[[i]] <- learn.VERSO.phylogenetic.tree(D = D, 
                                                            alpha = alpha[i], 
                                                            beta = beta[i], 
                                                            initialization = initialization, 
                                                            num_rs = num_rs, 
                                                            num_iter = num_iter, 
                                                            n_try_bs = n_try_bs, 
                                                            verbose = verbose)
            
        }

    }
    else {

        # parallel computation
        res_clusterEvalQ <- clusterEvalQ(parallel,library("Rfast"))
        clusterExport(parallel,varlist=c("D","alpha","beta","initialization","num_rs","num_iter","n_try_bs","verbose"),envir=environment())
        clusterExport(parallel,c("learn.VERSO.phylogenetic.tree","initialize.B","move.B","compute.C"),envir=environment())
        clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
        inference <- parLapply(parallel,seq_len(length(alpha)),function(x) {
            
            if(verbose) {
                message(paste('Performing inference for alpha =',paste0(alpha[x],collapse=" | "),'and beta =',paste0(beta[x],collapse=" | ")))
            }
            
            inference <- learn.VERSO.phylogenetic.tree(D = D, 
                                                       alpha = alpha[x], 
                                                       beta = beta[x], 
                                                       initialization = initialization, 
                                                       num_rs = num_rs, 
                                                       num_iter = num_iter, 
                                                       n_try_bs = n_try_bs, 
                                                       verbose = verbose)

        })
        
    }
    
    # close parallel
    if(close_parallel) {
        stopCluster(parallel)
    }

    # return the solution at maximum log-likelihood among the inferrend ones
    lik <- NULL
    for(i in seq_len(length(inference))) {
        lik <- c(lik,inference[[i]][["log_likelihood"]])
    }
    best <- which(lik==max(lik))[1]
    inference <- inference[[best]]
    B <- inference$B
    storage.mode(B) <- "integer"
    C <- inference$C
    log_likelihood <- inference[["log_likelihood"]]
    error_rates <- list(alpha=alpha[best],beta=beta[best])

    # compute corrected genotypes
    C_matrix <- array(0,c(nrow(D),ncol(B)))
    rownames(C_matrix) <- rownames(D)
    colnames(C_matrix) <- colnames(B)
    for(i in seq_len(nrow(C_matrix))) {
        C_matrix[i,which(rownames(B)==C[rownames(C_matrix)[i],"Genotype"])] <- 1
    }
    corrected_genotypes <- C_matrix %*% B
    corrected_genotypes <- corrected_genotypes[,colnames(D)]
    storage.mode(corrected_genotypes) <- "integer"

    # compute genotypes prevalence
    genotypes_prevalence <- array(NA,c(nrow(B),2))
    rownames(genotypes_prevalence) <- rownames(B)
    colnames(genotypes_prevalence) <- c("#Patients","Prevalence")
    patients <- table(C[,"Genotype"])
    prevalence <- patients/nrow(C)
    genotypes_prevalence[names(patients),"#Patients"] <- as.numeric(patients)
    genotypes_prevalence[names(prevalence),"Prevalence"] <- as.numeric(prevalence)

    # compute genotypes summary
    genotypes_summary <- list()
    i = 1
    for(genotypes in rownames(B)) {
        mut_list <- colnames(B)[B[i,]==1]
        genotypes_summary[[genotypes]] <- mut_list
        i = i + 1
    }

    # finally build VERSO phylogenetic tree
    phylogenetic_tree <- get.phylo(adjacency_matrix=as.adj.matrix(B),valid_genotypes=B,samples_attachments=C)

    # return results of the inference
    results <- list(B=B,C=C,phylogenetic_tree=phylogenetic_tree,corrected_genotypes=corrected_genotypes,genotypes_prevalence=genotypes_prevalence,genotypes_summary=genotypes_summary,log_likelihood=log_likelihood,error_rates=error_rates)
    return(results)

}
