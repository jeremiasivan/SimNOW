# functions for codes/1_sequence_simulation

# function: open msfile and convert locus trees into data.table
# required library: ape, data.table
f_extract_locus_tree <- function(fn_msfile, len_aln) {
    # open msfile
    ms <- data.table::fread(fn_msfile, col.names="cmd", sep="")
  
    # extract ms trees
    mstrees <- ms[grepl("^\\[", cmd)]
    if (nrow(mstrees) == 0) {
        mstrees <- ms[grepl("^\\(.+\\);$", cmd)]
        mstrees[1, cmd := paste0("[", len_aln, "]", mstrees[1, cmd])]
    }
    
    # convert ms trees into data.table
    mstrees[, c("empty","length","tree") := data.table::tstrsplit(cmd, '\\[|\\]', fixed = FALSE)]
    mstrees[, c("cmd", "empty") := NULL]
    mstrees$length <- as.numeric(mstrees$length)
    
    # extract topology and boundary of the locus trees
    n <- nrow(mstrees)
    start <- 1
    stop <- mstrees$length[1]
    
    # store the topology of the first locus
    first_tree <- ape::read.tree(text=mstrees$tree[1])
    first_tree$edge.length <- NULL
    topology <- ape::write.tree(first_tree)
        
    # store the topology of subsequent loci
    if (n > 1) {
        for (i in 2:n){
            start <- c(start,stop[length(stop)]+1)
            stop <- c(stop,stop[length(stop)]+mstrees$length[i])

            temp_tree <- ape::read.tree(text=mstrees$tree[i])
            temp_tree$edge.length <- NULL
            topology <- c(topology, ape::write.tree(temp_tree))
        }
    }

    # add new columns to the data.table
    mstrees[, c("start", "stop", "topology") := list(start, stop, topology)]
    setcolorder(mstrees, c("tree","length","start","stop","topology"))

    # match the naming convention of topologies
    uq_tops <- unique(mstrees$topology)
    ls_unroot <- NULL
    ls_uqtops <- NULL

    for (i in uq_tops) {
        is_equal <- FALSE
        
        for (j in ls_unroot) {
            # store the index of duplicated topology by referring to the existing entry
            if (ape::all.equal.phylo(ape::unroot(ape::read.tree(text = i)), ape::unroot(ape::read.tree(text = j)))) {
                is_equal <- TRUE
                ls_uqtops[[i]] <- j
                break
            }
        }
        
        # store the index of unique topology as a new entry
        if (!is_equal) {
            ls_unroot <- c(ls_unroot, i)
            ls_uqtops[[i]] <- i
        }
    }
    
    # update the topology based on the entries
    lookup <- match(mstrees$topology, names(ls_uqtops))
    mstrees[, "topology" := as.character(ls_uqtops[lookup])]

    return(mstrees)
}

# function: open mstrees and combine contiguous loci with the same topology
# required library: dplyr, data.table
f_extract_locus_topology <- function(mstrees) {
    # combine contiguous topologies together 
    t_summary <- mstrees %>%
        group_by(topology, group_run = data.table::rleid(topology)) %>%
        summarise(length=sum(length)) %>%
        arrange(group_run)
    t_summary$group_run <- NULL
    t_summary <- as.data.table(t_summary)
    
    # store the boundary of each topology 
    nt <- nrow(t_summary)
    start <- 1
    stop <- t_summary$length[1]

    if (nt > 1) {
        for (i in 2:nt){
            start <- c(start,stop[length(stop)]+1)
            stop <- c(stop,stop[length(stop)]+t_summary$length[i])
        }
    } 

    # add new columns to the data.table
    t_summary[, c("start", "stop") := list(start, stop)]

    return(t_summary)
}

# function: copy gap positions from source to target alignments
# required library: seqinr
f_copy_gaps <- function(fn_target_aln, fn_source_aln) {
    # read source and simulated alignments
    srcseq <- seqinr::read.fasta(fn_source_aln, whole.header=T)
    simseq <- seqinr::read.fasta(fn_target_aln, whole.header=T)
    
    # trim whitespaces from FASTA headers
    names(srcseq) <- sapply(names(srcseq), function(x){trimws(x, which="both")})
    names(simseq) <- sapply(names(simseq), function(x){trimws(x, which="both")})
    
    # extract FASTA headers
    srctaxa <- sort(names(srcseq))
    simtaxa <- sort(names(simseq))

    # return error message if the taxa do not match
    if (!all(srctaxa == simtaxa)) {
      err_msg <- "Error: source alignment has different taxa than the simulated one. Exited."
      return(list(seq=NULL, err=err_msg))
    }
    
    # loop over taxon and paste gaps
    for (t in srctaxa) {
      gaps <- grep("-", srcseq[[t]])
      simseq[[t]][gaps] <- "-"
    }

    # return target sequence with gaps
    return(list(seq=simseq, err=NULL))
}