# Perform a robust estimation by bootstrap of alpha coefficients until a given level of cosine similarity is reached.
"signaturesSignificance" <- function( x, beta, cosine_thr = 0.95, min_contribution = 0.05, pvalue_thr = 0.05, sparsify = FALSE, nboot = 100, num_processes = Inf, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)

    # check input parameters
    x <- as.matrix(x)
    beta <- as.matrix(beta)
    if(nboot<5) {
        warning("The minimum value of nboot can be 5...")
        nboot <- 5
    }

    # setting up parallel execution
    parallel <- NULL
    close_parallel <- FALSE
    if(is.na(num_processes) || is.null(num_processes)) {
        num_processes <- 1
    }
    else if(num_processes==Inf) {
        cores <- as.integer((detectCores()-1))
        num_processes <- min(cores,nboot)
    }
    else {
        num_processes <- min(num_processes,nboot)
    }

    if(verbose) {
        cat("Estimating the contribution of each signature to the fit with a total of",nboot,"bootstrap iterations...","\n")
        if(num_processes>1) {
            cat("Executing",num_processes,"processes in parallel...","\n")
        }
    }

    # performing inference
    if(num_processes==1) {

        # performing inference sequentially
        alpha <- list()
        for(boot_iteration in 1:nboot) {

            if(verbose) {
                cat(paste0("Performing iteration ",boot_iteration," out of ",nboot,"..."),"\n")
            }

            curr_alpha <- lapply(X=1:nrow(x),FUN=function(counts) {
                curr_patient <- x[counts,,drop=FALSE]
                curr_num_counts <- round(sum(curr_patient)*100)
                curr_counts_distribution <- as.numeric(curr_patient/rowSums(curr_patient))
                curr_boot_sampling <- table(sample(x=colnames(curr_patient),size=curr_num_counts,replace=TRUE,prob=curr_counts_distribution))
                curr_patient[1,] <- 0
                curr_patient[1,names(curr_boot_sampling)] <- as.numeric(curr_boot_sampling)/100
                return(sigs.significance(x=curr_patient,beta=beta,cosine_thr=cosine_thr,sparsify=sparsify))
            })
            alpha[[boot_iteration]] <- Reduce("rbind",curr_alpha)

        }

    }
    else {

        # performing inference in parallel
        parallel <- makeCluster(num_processes,outfile="")
        close_parallel <- TRUE
        res_clusterEvalQ <- clusterEvalQ(parallel,library("nnls",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE))
        res_clusterEvalQ <- clusterEvalQ(parallel,library("glmnet",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE))
        res_clusterEvalQ <- clusterEvalQ(parallel,library("lsa",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE))
        clusterExport(parallel,varlist=c("sigs.significance"),envir=environment())
        clusterExport(parallel,varlist=c("verbose","nboot"),envir=environment())
        clusterExport(parallel,varlist=c("x","beta","cosine_thr","sparsify"),envir=environment())
        clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
        rm(res_clusterEvalQ)
        gc(verbose=FALSE)
        alpha <- parLapply(parallel,1:nboot,function(boot_iteration) {

            if(verbose) {
                cat(paste0("Performing iteration ",boot_iteration," out of ",nboot,"..."),"\n")
            }

            curr_alpha <- lapply(X=1:nrow(x),FUN=function(counts) {
                curr_patient <- x[counts,,drop=FALSE]
                curr_num_counts <- round(sum(curr_patient)*100)
                curr_counts_distribution <- as.numeric(curr_patient/rowSums(curr_patient))
                curr_boot_sampling <- table(sample(x=colnames(curr_patient),size=curr_num_counts,replace=TRUE,prob=curr_counts_distribution))
                curr_patient[1,] <- 0
                curr_patient[1,names(curr_boot_sampling)] <- as.numeric(curr_boot_sampling)/100
                return(sigs.significance(x=curr_patient,beta=beta,cosine_thr=cosine_thr,sparsify=sparsify))
            })
            curr_alpha <- Reduce("rbind",curr_alpha)

            return(curr_alpha)

        })

    }

    # process results
    results <- alpha[[1]]
    if(length(alpha)>1) {
        for(i in 2:length(alpha)) {
            results <- results + alpha[[i]]
        }
        results <- results/length(alpha)
    }
    alpha_distibution <- alpha
    alpha <- results

    if(verbose) {
        cat("Estimating level of significance for each signature...","\n")
    }
    pvalues <- alpha
    pvalues[which(pvalues==0)] <- NA
    for(i in 1:nrow(pvalues)) {
        for(j in 1:ncol(pvalues)) {
            if(!is.na(pvalues[i,j])) {
                pvalues[i,j] <- NA
                curr_values = NULL
                for(k in 1:length(alpha_distibution)) {
                    curr_values <- c(curr_values,(alpha_distibution[[k]][i,]/sum(alpha_distibution[[k]][i,]))[j])
                }
                pvalues[i,j] <- wilcox.test(as.numeric(curr_values),alternative="greater",mu=min_contribution)$p.value
            }
        }
    }
    goodness_fit <- NULL
    for(i in 1:nrow(x)) {
        # estimate goodness of fit
        goodness_fit <- c(goodness_fit,as.numeric(cosine(as.numeric(x[i,]),as.numeric((alpha[i,]%*%beta)))))
    }
    names(goodness_fit) <- rownames(x)
    bootstrap <- list(estimate=alpha,pvalues=pvalues,goodness_fit=goodness_fit)

    if(verbose) {
        cat("Performing fit of alpha considering only signatures with significant contribution...","\n")
    }

    # perform final fit of alpha using only significant signatures
    alpha <- array(0,c(nrow(x),nrow(beta)))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    goodness_fit <- NULL
    for(i in 1:nrow(x)) {
        # perform fit
        curr_sigs <- names(which(bootstrap$pvalues[i,]<pvalue_thr))
        curr_alpha <- array(NA,c(1,length(curr_sigs)))
        curr_beta <- beta[curr_sigs,,drop=FALSE]
        if(sparsify==TRUE&&nrow(curr_beta)>1) {
            res <- cv.glmnet(t(curr_beta),as.vector(x[i,]),type.measure="mse",nfolds=10,nlambda=10,family="gaussian",lower.limits=0.00)
            curr_alpha[1,] <- as.numeric(res$glmnet.fit$beta[,ncol(res$glmnet.fit$beta)])
        }
        else {
            curr_alpha[1,] <- nnls(t(curr_beta),as.vector(x[i,]))$x
        }
        alpha[i,rownames(curr_beta)] <- as.numeric(curr_alpha)
        # estimate goodness of fit
        goodness_fit <- c(goodness_fit,as.numeric(cosine(as.numeric(x[i,]),as.numeric((alpha[i,]%*%beta)))))
    }
    names(goodness_fit) <- rownames(x)
    
    # close parallel
    if(close_parallel) {
        stopCluster(parallel)
    }
    rm(close_parallel)
    rm(parallel)
    gc(verbose=FALSE)

    # save results
    results <- list(alpha=alpha,beta=beta,goodness_fit=goodness_fit,bootstrap_estimates=bootstrap)

    # return the estimeted signatures contributions
    return(results)
    
}

# iteratively estimate alpha coefficients until a given level of cosine similarity is reached
"sigs.significance" <- function( x, beta, cosine_thr = 0.95, sparsify = FALSE ) {

    # initialization
    alpha <- array(NA,c(1,nrow(beta)))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    if(sparsify==TRUE&&nrow(beta)>1) {
        res <- cv.glmnet(t(beta),as.vector(x[1,]),type.measure="mse",nfolds=10,nlambda=10,family="gaussian",lower.limits=0.00)
        alpha[1,] <- as.numeric(res$glmnet.fit$beta[,ncol(res$glmnet.fit$beta)])
    }
    else {
        alpha[1,] <- nnls(t(beta),as.vector(x[1,]))$x
    }

    # iteratively include signatures into the fit until a given level of cosine similarity is reached
    sigs <- colnames(alpha[,which(alpha[1,]>0),drop=FALSE])[sort.int(alpha[,which(alpha[1,]>0)],decreasing=TRUE,index.return=TRUE)$ix]
    alpha <- array(0,c(1,nrow(beta)))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    for(i in 1:length(sigs)) {
        # consider current set of signatures and perform fit of alpha
        curr_alpha <- array(NA,c(1,i))
        curr_beta <- beta[sigs[1:i],,drop=FALSE]
        if(sparsify==TRUE&&nrow(curr_beta)>1) {
            res <- cv.glmnet(t(curr_beta),as.vector(x[1,]),type.measure="mse",nfolds=10,nlambda=10,family="gaussian",lower.limits=0.00)
            curr_alpha[1,] <- as.numeric(res$glmnet.fit$beta[,ncol(res$glmnet.fit$beta)])
        }
        else {
            curr_alpha[1,] <- nnls(t(curr_beta),as.vector(x[1,]))$x
        }
        # estimate goodness of fit
        curr_predicted_counts <- curr_alpha %*% curr_beta
        curr_goodness_fit <- as.numeric(cosine(as.numeric(x),as.numeric(curr_predicted_counts)))
        if(curr_goodness_fit>cosine_thr) {
            break;
        }
    }
    alpha[,rownames(curr_beta)] <- curr_alpha

    # return the estimated alpha coefficients
    return(alpha)

}
