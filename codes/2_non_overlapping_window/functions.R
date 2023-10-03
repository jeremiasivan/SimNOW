# functions for codes/2_non_overlapping_window

# function: generate window alignments
# required library: seqinr
f_window_alignment <- function(fasta_aln, window_size, len_window, dir_output) {
    start <- 1
    wi <- ceiling(log(len_window) / log(10)) + 1
    
    for (i in 1:len_window) {
        subfasta <- lapply(fasta_aln, function(x) x[seq(from = start, to = as.numeric(i*window_size))])
        df <- do.call(rbind,subfasta)
        
        subfasta <- setNames(split(df, seq(nrow(df))), rownames(df))
        seqinr::write.fasta(sequences=subfasta, names=names(subfasta),
                    file.out=paste0(dir_output,"window_",formatC(i,width=wi,format="d",flag="0"),".fa"),
                    nbchar = 100)
        
        start <- start + window_size
    }
}

# function: generate window trees
f_window_tree <- function(dir_aln, prefix, thread, outgroup, set_model, set_blmin, dna_model, window_size, dir_iqtree2) {
    iqtree_cmd <- paste(dir_iqtree2,
                        "-S", dir_aln,
                        "-pre", prefix,
                        "-T", thread,
                        "-cptime 1000000 --quiet -redo")
    
    # set outgroup
    if (!is.null(outgroup) && outgroup == ""){
        iqtree_cmd <- paste(iqtree_cmd, "-o", outgroup)
    }
    
    # set substitution model
    if (set_model){
        iqtree_cmd <- paste(iqtree_cmd, "-m", dna_model)
    }
    
    # set minimum branch length (i.e., minimum one substitution per window)
    if (set_blmin) {
        iqtree_cmd <- paste(iqtree_cmd, "-blmin", 1/window_size)
    }
    
    system(iqtree_cmd)
}

# function: annotate window trees
# required library: data.table, dplyr
f_annotate_window_tree <- function(fn_treefile, list_window, window_size, prefix, dir_summary){
    # output files
    fn_atsum <- paste0(dir_summary, prefix, ".atsum")
    fn_cnsum <- paste0(dir_summary, prefix, ".cnsum")
    fn_tops <- paste0(dir_summary, prefix, ".topsum")

    # open the treefile
    iqtree_file <- file(description = fn_treefile, open="r", blocking = TRUE)
  
    # initiate variables
    seq_len <- c()
    start <- 0
    
    # iterate through each window tree
    while (length(trl <- readLines(iqtree_file, n=1))) {
        tl <- ape::read.tree(text=trl)
        tl$edge.length <- NULL
        tree <- ape::write.tree(tl)
        
        # store the values in dataframe
        pl <- list_window[start/window_size+1]
        seq_len <- rbind(seq_len, c(pl,start+1,start+window_size,window_size,tree))
        
        start <- start+window_size
    }
    close(iqtree_file)
    
    # convert vector into data.table
    seq_len <- data.table::as.data.table(seq_len)
    data.table::setnames(seq_len, c("source", "start", "stop", "length", "topology"))
    data.table::fwrite(seq_len, file=fn_atsum, quote=F, sep="\t")
    
    # extract top topologies
    top_tops <- seq_len %>%
        group_by(topology) %>%
        summarise(count = n()) %>%
        arrange(desc(count)) %>%
        mutate(cum.percentage = round(cumsum(count)/nrow(seq_len),3))
    data.table::fwrite(top_tops, file=fn_tops, quote=F, sep="\t")

    # count consecutive/contiguous windows per topology
    count_contiguous <- seq_len %>%
        group_by(topology,
                group_run = data.table::rleid(topology)) %>%
        summarise(count = n()) %>%
        arrange(group_run)
    count_contiguous$group_run <- NULL
    data.table::fwrite(count_contiguous, file=fn_cnsum, quote=F, sep="\t")
}

# function: compare topology weight between ms and NOW results
# required library: data.table, dplyr
f_compare_topology_weight <- function(msfile, smfile, fn_output){
    # count the topology weight of locus trees
    ms_tw <- msfile %>%
        group_by(topology) %>%
        summarise(sumlen = sum(length), percentage = sumlen/params$ms_l * 100)
    
    # count the topology weight from non-overlapping windows
    sm_tw <- smfile %>%
        group_by(topology) %>%
        summarise(sumlen = sum(length), percentage = sumlen/params$ms_l * 100)
    
    # add topology weight for each unique ms topology
    cmptw <- c()
    for (j in 1:nrow(ms_tw)) {
        top_idx <- match(ms_tw$topology[j], sm_tw$topology)
        if (is.na(top_idx)) {
            cmptw <- rbind(cmptw, c(ms_tw$topology[j], round(ms_tw$percentage[j], 4), "NA"))
        } else {
            cmptw <- rbind(cmptw, c(ms_tw$topology[j], round(ms_tw$percentage[j], 4), round(sm_tw$percentage[top_idx], 4)))
        }
    }
    
    # add topology weight for each unique NOW topology
    for (j in subset(sm_tw$topology, !sm_tw$topology %in% mstree_uq)) {
        top_idx <- match(j, sm_tw$topology)
        cmptw <- rbind(cmptw, c(j, "NA", round(sm_tw$percentage[top_idx], 4)))
    }
    
    # convert vector to data.table
    cmptw <- data.table::as.data.table(cmptw)
    data.table::setnames(cmptw, c("topology", "ms", "now"))
    data.table::fwrite(cmptw, file=fn_output, quote=FALSE, sep="\t")
}

# function: create site assignment matrix between ms and NOW results
# required library: data.table
f_create_site_acc_matrix <- function(msfile, smfile, mstree_uq, fn_output) {
    # count the number of unique topology
    mstree_uq_ln <- length(mstree_uq)

    # create empty dataframe
    output <- data.table::data.table(c(mstree_uq,"NT"), matrix(0, nrow=mstree_uq_ln+1, ncol=mstree_uq_ln))
    data.table::setnames(output, c("topology", mstree_uq))
    
    # iterate through ms topology file
    for (k in 1:nrow(msfile)){
        # extract the topology based on ms topology switching
        start_idx <- tail(which(smfile$start <= msfile$start[k]), 1)
        stop_idx <- head(which(smfile$stop >= msfile$stop[k]), 1)
        
        # update lengths
        sites <- smfile[start_idx:stop_idx,]
        sites$start[1] <- msfile$start[k]
        sites$length[1] <- sites$stop[1] - sites$start[1] + 1
        
        last_idx <- nrow(sites)
        sites$stop[last_idx] <- msfile$stop[k]
        sites$length[last_idx] <- sites$stop[last_idx] - sites$start[last_idx] + 1
        
        # extract index of the topology
        topology_idx <- match(msfile$topology[k], colnames(output))
        
        # add number of window sites for the respective ms topology
        for (j in 1:nrow(sites)) {
            window_idx <- match(sites$topology[j], output$topology)
            
            if (is.na(window_idx)) {
                tempval <- as.numeric(output[mstree_uq_ln+1, ..topology_idx] + sites$length[j])
                output[mstree_uq_ln+1, eval(msfile$topology[k]) := tempval]
            } else {
                tempval <- as.numeric(output[window_idx, ..topology_idx] + sites$length[j])
                output[window_idx, eval(msfile$topology[k]) := tempval]
            }
        }
    }
    
    data.table::fwrite(output, file=fn_output, sep = "\t", quote = F)
}

# function: extract accuracy and information criterion score from window trees
# required library: data.table
f_extract_acc_and_ic <- function(fn_cmp, fn_cmptw, fn_iqtree){
    # open the site assignment matrix file
    cmp <- data.table::fread(fn_cmp)
    cmp <- as.matrix(cmp[,-1])
    
    # calculate the site accuracy
    acc <- sum(as.numeric(diag(cmp))) / sum(as.numeric(cmp)) * 100
    acc <- round(acc, 3)
    
    # open the topology weight file
    cmptw <- data.table::fread(fn_cmptw)

    # set NULL value as zero
    na_ms <- which(is.na(cmptw$ms))
    na_now <- which(is.na(cmptw$now))
    
    if (length(na_ms) > 0) {
        cmptw$ms[na_ms] <- 0
    }
    
    if (length(na_now) > 0) {
        cmptw$now[na_now] <- 0
    }
    
    # extract the RMSE of the topology distribution
    rmse <- mean((cmptw$ms - cmptw$now) ^ 2) %>% sqrt()
    
    # extract information criterion scores
    aic <- gsub("^.* ", "", system(paste("grep '^Akaike information criterion'",fn_iqtree), intern = T))
    aicc <- gsub("^.* ", "", system(paste("grep '^Corrected Akaike information criterion'",fn_iqtree), intern = T))
    bic <- gsub("^.* ", "", system(paste("grep '^Bayesian information criterion'",fn_iqtree), intern = T))

    return(c(acc,aic,aicc,bic,rmse))
}