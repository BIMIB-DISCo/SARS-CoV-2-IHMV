# build a phylogenetic tree from a variants tree (adjacency_matrix) and corrected genotypes (theoretical_genotypes)
"VERSOphylo" <- function( adjacency_matrix, theoretical_genotypes ) {
    
    # compute valid genotypes given the adjacency matrix
    graph <- graph_from_adjacency_matrix(adjacency_matrix)
    roots <- names(which(colSums(adjacency_matrix)==0))
    valid_genotypes <- NULL
    genotypes_names <- NULL
    for(r in roots) {
        curr_genotype <- rep(0,nrow(adjacency_matrix))
        curr_genotype[which(rownames(adjacency_matrix)==r)] <- 1
        valid_genotypes <- rbind(valid_genotypes,curr_genotype)
        genotypes_names <- c(genotypes_names,r)
        curr <- all_simple_paths(graph=graph,from=r)
        if(length(curr)>0) {
            for(p in 1:length(curr)) {
                curr_path <- as.character(vertex_attr(graph,index=as.vector(curr[[p]]))$name)
                curr_path_pos <- NULL
                for(cpp in curr_path) {
                    curr_path_pos <- c(curr_path_pos,which(rownames(adjacency_matrix)==cpp))
                }
                curr_genotype <- rep(0,nrow(adjacency_matrix))
                curr_genotype[curr_path_pos] <- 1
                valid_genotypes <- rbind(valid_genotypes,curr_genotype)
                genotypes_names <- c(genotypes_names,cpp)
            }
        }
    }
    rownames(valid_genotypes) <- genotypes_names
    colnames(valid_genotypes) <- rownames(adjacency_matrix)

    # add Root/Reference genotype
    valid_genotypes <- valid_genotypes[colnames(valid_genotypes),]
    valid_genotypes <- rbind(rep(0,ncol(valid_genotypes)),valid_genotypes)
    rownames(valid_genotypes)[1] <- "Root"
    
    # compute Manhattan distance among valid genotypes
    distance_genotypes <- as.matrix(dist(valid_genotypes,method="manhattan"))

    # get samples attachments given theoretical and valid genotypes
    samples_attachment <- array(NA,c(nrow(theoretical_genotypes),1))
    rownames(samples_attachment) <- rownames(theoretical_genotypes)
    colnames(samples_attachment) <- "Genotype"
    for(gtG in rownames(valid_genotypes)) {
        curr_attachment <- which(apply(theoretical_genotypes,1,function(x) {
            all(as.numeric(x)==as.numeric(valid_genotypes[gtG,]))
        }))
        samples_attachment[names(curr_attachment),"Genotype"] <- gtG
    }

    # consider the variants tree (adjacency_matrix) to build inner nodes structure of the phylogenetic tree
    adjacency_matrix_with_reference <- rbind(rep(0,ncol(adjacency_matrix)),adjacency_matrix)
    adjacency_matrix_with_reference <- cbind(rep(0,nrow(adjacency_matrix_with_reference)),adjacency_matrix_with_reference)
    adjacency_matrix_with_reference[1,which(colSums(adjacency_matrix_with_reference)==0)[-1]] <- 1
    rownames(adjacency_matrix_with_reference)[1] <- "Root"
    colnames(adjacency_matrix_with_reference)[1] <- "Root"
    edges <- which(adjacency_matrix_with_reference==1,arr.ind=TRUE)
    parents_list <- rownames(adjacency_matrix_with_reference)[as.numeric(edges[,"row"])]
    children_list <- rownames(adjacency_matrix_with_reference)[as.numeric(edges[,"col"])]
    edges <- cbind(parents_list,children_list)
    leaf_nodes <- rownames(adjacency_matrix_with_reference)[as.numeric(which(rowSums(adjacency_matrix_with_reference)==0))]
    observed_attachments <- unique(samples_attachment[,"Genotype"])
    hidden_nodes <- leaf_nodes[which(!leaf_nodes%in%observed_attachments)]
    nodes_list <- unique(c(parents_list,children_list))
    hidden_nodes <- as.numeric(which(nodes_list%in%hidden_nodes))
    for(mli in 1:length(nodes_list)) {
        curr_mli <- which(parents_list==nodes_list[mli])
        if(length(curr_mli)>0) {
            parents_list[curr_mli] <- mli
        }
        curr_mli <- which(children_list==nodes_list[mli])
        if(length(curr_mli)>0) {
            children_list[curr_mli] <- mli
        }
    }
    edges_weights <- NULL
    for(egs in 1:nrow(edges)) {
        edges_weights <- c(edges_weights,distance_genotypes[edges[egs,1],edges[egs,2]])
    }
    Nnode <- length(nodes_list)
    edges <- cbind(as.numeric(parents_list),as.numeric(children_list))
    colnames(edges) <- c("Parent","Child")
    
    # now consider samples attachments
    attachments <- samples_attachment[,"Genotype"]
    tip.label <- names(attachments)
    tip.label <- c(tip.label,"Reference")
    if(length(hidden_nodes)>0) {
        tip.label <- c(tip.label,paste0("verso_hidden_genotype_",1:length(hidden_nodes)))
    }
    overhead <- length(tip.label)
    edges <- edges + overhead
    if(length(hidden_nodes)>0) {
        hidden_nodes <- hidden_nodes + overhead
    }
    for(att in 1:length(attachments)) {
        edges <- rbind(edges,c((which(nodes_list==as.character(attachments[att]))+overhead),att))
        edges_weights <- c(edges_weights,0)
    }

    # finally add a sample attached to Root representing Reference and any attachment to hidden genotypes
    edges <- rbind(edges,c((1+overhead),(length(attachments)+1)))
    edges_weights <- c(edges_weights,0)
    if(length(hidden_nodes)>0) {
        for(att in 1:length(hidden_nodes)) {
            edges <- rbind(edges,c(hidden_nodes[att],which(tip.label==paste0("verso_hidden_genotype_",att))))
            edges_weights <- c(edges_weights,0)
        }
    }
    rownames(edges) <- 1:nrow(edges)

    # build a phylogenetic tree from the results by VERSO
    phylogenetic_tree <- list(edge=edges,tip.label=tip.label,Nnode=Nnode,edge.length=edges_weights,node.label=nodes_list,root.edge=0)
    class(phylogenetic_tree) <- "phylo"
    
    # return VERSO phylogenetic tree
    return(phylogenetic_tree)

}

# convert B to an adjacency matrix
"as.adj.matrix" <- function( B ) {

    # create the data structure where to save the adjacency matrix obtained from B
    adj_matrix <- array(0L,dim(B))
    rownames(adj_matrix) <- colnames(B)
    colnames(adj_matrix) <- colnames(B)

    # set arcs in the adjacency matrix
    for(i in 1:(nrow(B)-1)) {
        for(j in ((i+1):nrow(B))) {
            if(all(B[i,1:i]==B[j,1:i])&&(sum(B[i,])==(sum(B[j,])-1))) {
                adj_matrix[i,j] <- 1
            }
        }
    }

    # return the adjacency matrix obtained from B
    return(adj_matrix)

}
