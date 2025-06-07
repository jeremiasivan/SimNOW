# function: create per-window aligments
# required package: seqinr
f_perwindow_run <- function(dir_perwindow, fasta, window_name, start, end) {
  subfasta <- lapply(fasta, function(x) x[seq(from = start, to = end)])
  subfasta <- do.call(rbind,subfasta)
  subfasta <- setNames(split(subfasta, seq(nrow(subfasta))), rownames(subfasta))

  # write FASTA file
  fn_out <- paste0(dir_perwindow, "/", window_name, ".fa")
  seqinr::write.fasta(sequences=subfasta, names=names(subfasta), file.out=fn_out, nbchar=100)
}

# function: create window tree with duplicate seqs
f_iqtree2_single <- function(input, outgroup, window_size, setblmin, setkeepident, setmodel, dna_model, bs_type, bs, dir_iqtree2) {
  iqtree_cmd <- paste(dir_iqtree2,
                      "-s", input,
                      "-T 1 --quiet -redo")
  
  if (!is.null(outgroup) && !outgroup == ""){
    iqtree_cmd <- paste(iqtree_cmd, "-o", outgroup)
  }
  
  if (setblmin) {
    iqtree_cmd <- paste(iqtree_cmd, "-blmin", 1/window_size)
  }
  
  if (setkeepident) {
    iqtree_cmd <- paste(iqtree_cmd, "-keep-ident")
  }

  if (setmodel) {
    iqtree_cmd <- paste(iqtree_cmd, "-m", dna_model)
  }

  if (!is.null(bs_type) && bs_type != "") {
    if (tolower(bs_type) == "ufboot") {
      iqtree_cmd <- paste(iqtree_cmd, "-bb", bs)
    } else if (tolower(bs_type) == "nonparametric") {
      iqtree_cmd <- paste(iqtree_cmd, "-b", bs)
    }
  }
  
  system(iqtree_cmd)
}

# function: extract window tree statistics
# required package: ape
f_window_tree_statistics <- function(fn_iqtree, fn_treefile, bootstrap_type, min_branch_support) {
  # extract log-likelihood and number of free parameters
  logl <- gsub(" \\(.*\\)$", "", system(paste("grep '^Log-likelihood of the tree'",fn_iqtree), intern = T))
  logl <- as.numeric(gsub("^.* ", "", logl))
  freeparams <- as.numeric(gsub("^.* ", "", system(paste("grep '^Number of free parameters'",fn_iqtree), intern = T)))
  
  # read treefile
  tre <- readLines(fn_treefile)

  tl <- ape::read.tree(text=tre)
  tl$edge.length <- NULL

  if (!is.null(bootstrap_type) && bootstrap_type != "") {
    bl <- subset(tl$node.label, tl$node.label != "")

    if (!is.null(bl) && mean(as.numeric(bl)) >= min_branch_support) {
      tl$node.label <- NULL
      return(list(logl=logl, freeparams=freeparams, tree=ape::write.tree(tl)))
    }

    return(list(logl=logl, freeparams=freeparams, tree=NULL))
  }

  return(list(logl=logl, freeparams=freeparams, tree=ape::write.tree(tl)))
}

# function: extract AIC from iqtree
f_window_tree_aic <- function(fn_iqtree) {
  # extract log-likelihood and number of free parameters
  logl <- gsub(" \\(.*\\)$", "", system(paste("grep '^Log-likelihood of the tree'",fn_iqtree), intern = T))
  logl <- as.numeric(gsub("^.* ", "", logl))
  freeparams <- as.numeric(gsub("^.* ", "", system(paste("grep '^Number of free parameters'",fn_iqtree), intern = T)))

  return(list(logl=logl, freeparams=freeparams))
}

# function: extract accuracy from window trees
# required library: data.table
f_extract_acc <- function(fn_cmp, fn_cmptw){
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

    return(list(acc=acc, rmse=rmse))
}
