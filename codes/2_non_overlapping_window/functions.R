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