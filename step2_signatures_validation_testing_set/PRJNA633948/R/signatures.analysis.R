# perform signatures discovery and rank estimation for a range of K signatures given a set of observations x
"signatures.decomposition" <- function( x, K, nmf_runs = 1000, num_processes = Inf, seed = NULL ) {
    
    # set the seed
    set.seed(seed)

    # setting up parallel execution
    pbackend <- NULL
    close_parallel <- FALSE
    if(is.na(num_processes) || is.null(num_processes)) {
        num_processes <- 1
    }
    else if(num_processes==Inf) {
        cores <- as.integer((detectCores()-1))
        num_processes <- min(cores,nmf_runs)
    }
    else {
        num_processes <- min(num_processes,nmf_runs)
    }
    if(num_processes==1) {
        pbackend <- "seq"
    }
    else {
        pbackend <- makeCluster(num_processes)
        res_clusterEvalQ <- clusterEvalQ(pbackend,library("NMF",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE))
        res_clusterEvalQ <- clusterEvalQ(pbackend,library("nnls",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE))
        clusterExport(pbackend,varlist=c("nmf.seed","nmf.nnls"),envir=environment())
        clusterSetRNGStream(pbackend,iseed=round(runif(1)*100000))
        close_parallel <- TRUE
        rm(res_clusterEvalQ)
        gc(verbose=FALSE)
    }

    cat("Performing signatures discovery and rank estimation...","\n")
    cat("Executing",num_processes,"processes in parallel...","\n")

    # perform signatures discovery and rank estimation
    alpha <- list()
    beta <- list()
    measures <- NULL
    for(i in 1:length(K)) {

        cat(paste0("Performing inference for K=",K[i],"..."),"\n")

        # perform the inference for current K
        results <- NULL
        while(is.null(results)) {
            results <- tryCatch({
                results <- nmf(x=x,rank=K[i],method=nmf.nnls,seed=nmf.seed,rng=round(runif(1)*10000),nrun=nmf_runs,.pbackend=pbackend)
                gc(verbose=FALSE)
                results
            }, error = function(e) {
                cat(paste0("An error has occurred: ",e$message),"\n")
                gc(verbose=FALSE)
                NULL
            }, finally = {
                gc(verbose=FALSE)
            })
        }

        alpha[[paste0(K[i],"_signatures")]] <- basis(results)
        beta[[paste0(K[i],"_signatures")]] <- coef(results)

        # compute and save quality measures
        curr_rss <- rss(results,x)
        curr_evar <- evar(results,x)
        curr_silhouette_alpha <- tryCatch({
            mean(silhouette(results,what="features")[,"sil_width"])
        }, error = function(e) {
            NA
        })
        curr_silhouette_beta <- tryCatch({
            mean(silhouette(results,what="samples")[,"sil_width"])
        }, error = function(e) {
            NA
        })
        curr_sparseness_alpha <- sparseness(as.vector(alpha[[paste0(K[i],"_signatures")]]))
        curr_sparseness_beta <- sparseness(as.vector(beta[[paste0(K[i],"_signatures")]]))
        if(nmf_runs>1) {
            curr_silhouette_consensus <- tryCatch({
                mean(silhouette(results,what="chc")[,"sil_width"])
            }, error = function(e) {
                NA
            })
            curr_measures <- matrix(c(K[i],curr_rss,curr_evar,curr_silhouette_alpha,curr_silhouette_beta,curr_sparseness_alpha,curr_silhouette_beta,cophcor(results),dispersion(results),curr_silhouette_consensus),nrow=1)
            colnames(curr_measures) <- c("Rank","Residual_Sum_Squares","Explained_Variance","Silhouette_Alpha","Silhouette_Beta","Sparseness_Alpha","Sparseness_Beta","Cophenetic_Coefficient","Dispersion_Coefficient","Silhouette_Consensus")
        }
        else {
            curr_measures <- matrix(c(K[i],curr_rss,curr_evar,curr_silhouette_alpha,curr_silhouette_beta,curr_sparseness_alpha,curr_silhouette_beta),nrow=1)
            colnames(curr_measures) <- c("Rank","Residual_Sum_Squares","Explained_Variance","Silhouette_Alpha","Silhouette_Beta","Sparseness_Alpha","Sparseness_Beta")
        }
        measures <- rbind(measures,curr_measures)

    }
    rownames(measures) <- 1:nrow(measures)

    # save results
    results <- list(alpha=alpha,beta=beta,measures=measures)
    
    # close pbackend
    if(close_parallel) {
        stopCluster(pbackend)
    }

    # return the discovered signatures
    return(results)
    
}

# perform the assessment of different signatures.decomposition solutions by cross validation given a set of observations x and discovered signatures beta
"signatures.cross.validation" <- function( x, beta, cross_validation_entries = 0.01, cross_validation_iterations = 5, cross_validation_repetitions = 1000, num_processes = Inf, seed = NULL ) {
    
    # set the seed
    set.seed(seed)

    # setting up parallel execution
    parallel <- NULL
    close_parallel <- FALSE
    if(is.na(num_processes) || is.null(num_processes)) {
        num_processes <- 1
    }
    else if(num_processes==Inf) {
        cores <- as.integer((detectCores()-1))
        num_processes <- min(cores,cross_validation_repetitions)
    }
    else {
        num_processes <- min(num_processes,cross_validation_repetitions)
    }

    cat("Estimating the optimal number of signatures with a total of",cross_validation_repetitions,"cross validation repetitions...","\n")
    cat("Executing",num_processes,"processes in parallel...","\n")

    # structure to save the results
    cv_estimates <- array(NA,c(cross_validation_repetitions,length(beta)))
    rownames(cv_estimates) <- paste0("Repetition_",1:cross_validation_repetitions)
    colnames(cv_estimates) <- 1:ncol(cv_estimates)

    # perform a total of cross_validation_repetitions repetitions of cross validation
    valid_entries <- which(x>0,arr.ind=TRUE)

    # performing inference
    if(num_processes==1) {

        for(cv_repetitions in 1:cross_validation_repetitions) {

            cat(paste0("Performing repetition ",cv_repetitions," out of ",cross_validation_repetitions,"..."),"\n")

            # randomly set the cross validation entries for the current iteration
            cv_entries <- valid_entries[sample(1:nrow(valid_entries),size=round(nrow(valid_entries)*cross_validation_entries),replace=FALSE),]
            
            # consider all the possible values of K
            cont <- 0
            for(num_signs in 1:length(beta)) {

                k <- nrow(beta[[num_signs]])
                if(cv_repetitions==1) {
                    colnames(cv_estimates)[num_signs] <- paste0(k,"_signatures")
                }

                cat(paste0("Performing estimation for K=",k,"..."),"\n")
                
                # repeat the estimation for a number of cross_validation_iterations
                x_cv <- x
                for(cv_iteration in 1:cross_validation_iterations) {

                    cat(paste0("Performing cross validation iteration ",cv_iteration," out of ",cross_validation_iterations,"..."),"\n")

                    # set a percentage of cross_validation_entries entries to 0 in order to perform cross validation
                    if(cv_iteration==1) {
                        x_cv[cv_entries] <- 0
                    }
                    else {
                        predicted_counts <- curr_results[["alpha"]]%*%curr_results[["beta"]]
                        x_cv[cv_entries] <- predicted_counts[cv_entries]
                    }

                    # columns in x_cv with all 0s are not allowed
                    if(any(colSums(x_cv)==0)) {
                        invalid_cols <- as.numeric(which(colSums(x_cv)==0))
                        for(inv_cols in invalid_cols) {
                            x_cv[sample(1:length(x_cv[,inv_cols]),size=1),inv_cols] <- 1e-05
                        }
                    }

                    # perform the inference
                    curr_results <- nmf.fit(x=x_cv,beta=beta[[num_signs]])

                }

                # save an estimate of the current solution
                curr_predicted_counts <- curr_results[["alpha"]]%*%curr_results[["beta"]]
                curr_true_considered_counts <- as.vector(x[cv_entries])
                curr_predicted_considered_counts <- as.vector(curr_predicted_counts[cv_entries])
                error <- mean((curr_true_considered_counts-curr_predicted_considered_counts)^2)
                cv_estimates[cv_repetitions,paste0(k,"_signatures")] <- error

                cont <- cont + 1
                cat("Progress",paste0(round((cont/length(beta))*100,digits=3),"%..."),"\n")

            }

        }

        # compute mean and median values of estimated cross validation error
        cv_mean <- NULL
        cv_median <- NULL
        for(i in 1:ncol(cv_estimates)) {
            cv_mean <- c(cv_mean,mean(cv_estimates[,i]))
            cv_median <- c(cv_median,median(cv_estimates[,i]))
        }
        cv_summary <- cbind(cv_mean,cv_median)
        rownames(cv_summary) <- colnames(cv_estimates)
        colnames(cv_summary) <- c("Mean cross-validation error","Median cross-validation error")

    }
    else {

        parallel <- makeCluster(num_processes,outfile="")
        close_parallel <- TRUE
        res_clusterEvalQ <- clusterEvalQ(parallel,library("nnls",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE))
        clusterExport(parallel,varlist=c("nmf.fit"),envir=environment())
        clusterExport(parallel,varlist=c("cross_validation_repetitions","cross_validation_entries"),envir=environment())
        clusterExport(parallel,varlist=c("cross_validation_iterations","valid_entries","beta","x"),envir=environment())
        clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
        curr_results <- parLapply(parallel,1:cross_validation_repetitions,function(cv_repetitions) {

            cat(paste0("Performing repetition ",cv_repetitions," out of ",cross_validation_repetitions,"..."),"\n")

            # randomly set the cross validation entries for the current iteration
            cv_entries <- valid_entries[sample(1:nrow(valid_entries),size=round(nrow(valid_entries)*cross_validation_entries),replace=FALSE),]
            
            # consider all the possible values of K
            cv_errors <- rep(NA,length(beta))
            cv_names <- NULL
            for(num_signs in 1:length(beta)) {

                if(cv_repetitions==1) {
                    cv_names <- c(cv_names,paste0(nrow(beta[[num_signs]]),"_signatures"))
                }

                # repeat the estimation for a number of cross_validation_iterations
                x_cv <- x
                for(cv_iteration in 1:cross_validation_iterations) {

                    # set a percentage of cross_validation_entries entries to 0 in order to perform cross validation
                    if(cv_iteration==1) {
                        x_cv[cv_entries] <- 0
                    }
                    else {
                        predicted_counts <- curr_results[["alpha"]]%*%curr_results[["beta"]]
                        x_cv[cv_entries] <- predicted_counts[cv_entries]
                    }

                    # columns in x_cv with all 0s are not allowed
                    if(any(colSums(x_cv)==0)) {
                        invalid_cols <- as.numeric(which(colSums(x_cv)==0))
                        for(inv_cols in invalid_cols) {
                            x_cv[sample(1:length(x_cv[,inv_cols]),size=1),inv_cols] <- 1e-05
                        }
                    }

                    # perform the inference
                    curr_results <- nmf.fit(x=x_cv,beta=beta[[num_signs]])

                }

                # save an estimate of the current solution
                curr_predicted_counts <- curr_results[["alpha"]]%*%curr_results[["beta"]]
                curr_true_considered_counts <- as.vector(x[cv_entries])
                curr_predicted_considered_counts <- as.vector(curr_predicted_counts[cv_entries])
                error <- mean((curr_true_considered_counts-curr_predicted_considered_counts)^2)
                cv_errors[num_signs] <- error

            }
            names(cv_errors) <- cv_names

            return(cv_errors)

        })

        # save the results from parallel computation
        colnames(cv_estimates) <- names(curr_results[[1]])
        for(par_res in 1:length(curr_results)) {
            cv_estimates[par_res,] <- curr_results[[par_res]]
        }

        # compute mean and median values of estimated cross validation error
        cv_mean <- NULL
        cv_median <- NULL
        for(i in 1:ncol(cv_estimates)) {
            cv_mean <- c(cv_mean,mean(cv_estimates[,i]))
            cv_median <- c(cv_median,median(cv_estimates[,i]))
        }
        cv_summary <- cbind(cv_mean,cv_median)
        rownames(cv_summary) <- colnames(cv_estimates)
        colnames(cv_summary) <- c("Mean cross-validation error","Median cross-validation error")

    }

    # save results
    results <- list(estimates=cv_estimates,summary=cv_summary)
    
    # close parallel
    if(close_parallel) {
        stopCluster(parallel)
    }

    # return results of cross validation
    return(results)
    
}

# initialize alpha and beta for nmf.nnls function
"nmf.seed" <- function( model, target ) {

    # initialize alpha with an empty matrix
    alpha <- array(NA,c(nrow(target),nbasis(model)))
    rownames(alpha) <- rownames(target)
    colnames(alpha) <- paste0("S",1:ncol(alpha))

    # randomly initialize beta
    beta <- matrix(runif(nbasis(model)*6),nrow=nbasis(model),ncol=6)
    beta <- (beta/rowSums(beta)) # beta rows (i.e., signatures) must sum to 1
    if(any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))),] <- 0
    }
    rownames(beta) <- colnames(alpha)
    colnames(beta) <- colnames(target)

    # update results
    basis(model) <- alpha
    coef(model) <- beta

    # return updated model
    return(model)

}

# perform NMF by Non-negative least squares
"nmf.nnls" <- function( x, seed ) {

    # initialization
    alpha <- basis(seed) # exposures matrix
    beta <- coef(seed) # signatures matrix
    n <- nrow(x) # n is the number of observations in x, i.e., the samples
    J <- ncol(x) # J is the number of categories, i.e., the 6 contexts
    K <- nrow(beta) # K is the number of signatures to be fitted

    # iteratively fit alpha and beta by Non-negative least squares (nnls)
    for(i in 1:20) {

        # update alpha, beta is kept fixed
        for(j in 1:n) {
            alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
        }

        # update beta, alpha is kept fixed
        for(k in 1:J) {
            beta[,k] <- nnls(alpha,as.vector(x[,k]))$x
        }

    }

    # normalize beta and perform final fit of alpha
    beta <- (beta/rowSums(beta))
    if(any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))),] <- 0
    }
    for(j in 1:n) {
        alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
    }
    colnames(alpha) <- paste0("S",1:ncol(alpha))
    rownames(beta) <- colnames(alpha)

    # update results
    basis(seed) <- alpha
    coef(seed) <- beta

    # return updated seed
    return(seed)

}

# perform fit of a given NMF solution by Non-negative least squares
"nmf.fit" <- function( x, beta ) {

    # initialization
    alpha <- array(NA,c(nrow(x),nrow(beta)))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    beta <- (beta/rowSums(beta))
    if(any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))),] <- 0
    }
    n <- nrow(x) # n is the number of observations in x, i.e., the samples
    J <- ncol(x) # J is the number of categories, i.e., the 6 contexts
    
    # iteratively fit alpha and beta by Non-negative least squares (nnls)
    for(i in 1:20) {

        # update alpha, beta is kept fixed
        for(j in 1:n) {
            alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
        }

        # update beta, alpha is kept fixed
        for(k in 1:J) {
            beta[,k] <- nnls(alpha,as.vector(x[,k]))$x
        }

    }

    # normalize beta and perform final fit of alpha
    beta <- (beta/rowSums(beta))
    if(any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))),] <- 0
    }
    for(j in 1:n) {
        alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
    }

    # return the results
    results <- list(alpha=alpha,beta=beta)
    return(results)

}
