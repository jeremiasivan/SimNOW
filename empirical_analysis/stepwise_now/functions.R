# functions for codes_empirical/steowise_now

# function: calculate FASTA alignment length
f_get_alignment_length <- function(filepath) {
  # set up variables
  alignment_length <- 0
  first_sequence <- FALSE
  
  # open the FASTA file
  con <- file(filepath, "r")
  
  # read the file line by line
  while (length(line <- readLines(con, n = 1)) > 0) {
    # check if the header is the first in the alignment
    if (grepl("^>+", line)) {
      if (first_sequence) {
        break
      }

      # update variable
      first_sequence <- TRUE

    } else {
      # update the alignment length
      alignment_length <- alignment_length + nchar(line)
    }
  }
  
  # close the file connection
  close(con)
  
  return(alignment_length)
}

# function: create data.frame for all windows
f_create_perwindow_df <- function(numw, fasta_len, window_size) {
  iterator <- 1:numw
  width <- ceiling(log(numw) / log(10)) + 1

  # set up start and end positions
  start_pos <- seq(1, fasta_len - window_size + 1, by = window_size)
  end_pos <- start_pos + window_size - 1

  # create a data.frame
  df_windows <- data.frame(window_name = paste0("window_", formatC(iterator, width=width, format="d", flag="0")),
                          start = start_pos, end = end_pos)
  return(df_windows)
}

# function: create per-window aligments
# required package: seqinr, doSNOW
f_perwindow_run <- function(dir_perwindow, fasta, numw, fasta_len, window_size, thread) {
  # generate data.frame for all windows
  df_windows <- f_create_perwindow_df(numw, fasta_len, window_size)

  # remove all alignments in the folder
  unlink(paste0(dir_perwindow,"*.fa"))
  
  # create doSNOW cluster
  nwcl <- makeCluster(thread)
  doSNOW::registerDoSNOW(nwcl)

  # generate window alignments
  foreach (i = 1:nrow(df_windows)) %dopar% {
    subfasta <- lapply(fasta, function(x) x[seq(from = df_windows$start[i], to = df_windows$end[i])])
    subfasta <- do.call(rbind,subfasta)
    subfasta <- setNames(split(subfasta, seq(nrow(subfasta))), rownames(subfasta))

    # write FASTA file
    fn_out <- paste0(dir_perwindow, "/", df_windows$window_name[i], ".fa")
    seqinr::write.fasta(sequences=subfasta, names=names(subfasta), file.out=fn_out, nbchar=100)
  }

  stopCluster(nwcl)
}

# function: identify if window alignment is informative or not
f_filter_uninformative_window <- function(df_seq, len_taxa, min_informative_sites) {
  df_seq <- as.data.frame(df_seq)

  # first filter: remove alignment with incomplete taxa
  is_any_all_gaps <- any(apply(df_seq, 1, function(x) all(x == "-")))
  if (nrow(df_seq) != len_taxa || is_any_all_gaps) {
    return(list(is_informative=FALSE, err_msg="Error: one or more taxa has all gaps"))
  }

  # second filter: remove aligment with all constant sites
  is_all_constant <- all(apply(df_seq, 2, function(x) all(x == x[1])))
  if (is_all_constant) {
    return(list(is_informative=FALSE, err_msg="Error: alignment consists of all constant sites"))
  }

  # # third filter: unique sequences
  # n_unique_seq <- nrow(unique(df_seq))
  # if (n_unique_seq != len_taxa) {
  #   return(list(is_informative=FALSE, err_msg="Error: some alignments are duplicate"))
  # }

  # fourth filter: minimum informative sites
  count_constant <- sum(apply(df_seq, 2, function(x) all(x == x[1])))
  n_informative_sites <- ncol(df_seq) - count_constant
  if (n_informative_sites < min_informative_sites) {
    return(list(is_informative=FALSE, err_msg=paste("Error: alignment only has", n_informative_sites, "non-constant sites")))
  }

  return(list(is_informative=TRUE, err_msg=""))
}

# function: create window alignment
# required packages: data.table, doSNOW, dplyr, seqinr
f_generate_perwindowsum <- function(dir_perwindow, fasta_len, window_size, len_taxa, min_informative_sites, thread) {
  # generate data.frame for all windows
  numw <- fasta_len/window_size
  df_windows <- f_create_perwindow_df(numw, fasta_len, window_size)

  # create doSNOW cluster
  nwcl <- makeCluster(thread)
  doSNOW::registerDoSNOW(nwcl)

  # export function
  snow::clusterExport(nwcl, "f_filter_uninformative_window")

  # generate window alignments
  ls_windows <- foreach (i = 1:nrow(df_windows), .combine='c') %dopar% {
    # set up variables
    file_fasta <- paste0(dir_perwindow, "/", df_windows$window_name[i], ".fa")

    # convert sequence to data.frame
    fasta <- seqinr::read.fasta(file_fasta, whole.header=T)
    fasta <- as.data.frame(do.call(rbind, fasta))
    output <- f_filter_uninformative_window(fasta, len_taxa, min_informative_sites)

    return(list(list(name=df_windows$window_name[i], start=df_windows$start[i], end=df_windows$end[i], is_informative=output$is_informative, err_msg=output$err_msg)))
  }

  stopCluster(nwcl)

  # convert to data.table
  df_output <- data.table::as.data.table(do.call(rbind, ls_windows), fill=TRUE) %>% arrange(name)

  return(df_output)
}

# function: create window tree
f_iqtree2_single <- function(input, outgroup, setblmin, setmodel, dna_model, bs_type, bs, dir_iqtree2) {
  iqtree_cmd <- paste(dir_iqtree2,
                      "-s", input,
                      "-T 1 --quiet -redo")
  
  if (!is.null(outgroup) && !outgroup == ""){
    iqtree_cmd <- paste(iqtree_cmd, "-o", outgroup)
  }
  
  if (setblmin) {
    iqtree_cmd <- paste(iqtree_cmd, "-blmin", 1/i)
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
  }

  return(list(logl=logl, freeparams=freeparams, tree=NULL))
}

# function: summary across window trees
# required package: data.table, dplyr
f_window_trees_summary <- function(ls_statistics, fn_uqtops) {
  # extract AIC
  logl <- sum(sapply(ls_statistics, function(x) x[[1]]))
  freeparams <- sum(sapply(ls_statistics, function(x) x[[2]]))
  aic <- (2 * freeparams) - (2 * logl)

  tree <- unlist(sapply(ls_statistics, function(x) {x[[3]]}))

  # calculate the topology distribution
  df_topology <- data.table::as.data.table(tree)
  data.table::setnames(df_topology, "topology")
  sorted_topology <- df_topology %>%
    group_by(topology) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) %>%
    mutate(cum.percentage=round(cumsum(n)/sum(n)*100,3))
  data.table::fwrite(sorted_topology, file=fn_uqtops, quote=F, sep="\t")

  return(list(aic=aic, tree=tree))
}

# functions: calculate AIC from iqtree file
f_calculate_aic_from_iqtree <- function(ls_iqtree) {
  total_logl <- 0
  total_freeparams <- 0

  for (fn_iqtree in ls_iqtree) {
    # extract log-likelihood and number of free parameters
    logl <- gsub(" \\(.*\\)$", "", system(paste("grep '^Log-likelihood of the tree'",fn_iqtree), intern = T))
    logl <- as.numeric(gsub("^.* ", "", logl))
    freeparams <- as.numeric(gsub("^.* ", "", system(paste("grep '^Number of free parameters'",fn_iqtree), intern = T)))

    # update variables
    total_logl <- total_logl + logl
    total_freeparams <- total_freeparams + freeparams
  }

  aic <- (2 * total_freeparams) - (2 * total_logl)
  return(aic)
}

# functions: plot delta AIC across chromosome
f_chromosomal_delta_aic <- function(df_aic_sum, long_wsize, short_wsize) {
  # add two more columns for visualization
  df_delta_aic_plot <- df_aic_sum %>%
      mutate(mid=floor((start+stop)/2), bg_color=ifelse(delta_aic=="NA","red","white"))
  
  # visualization
  plot <- ggplot(df_delta_aic_plot, aes(x=mid, y=delta_aic, xmin=min(start), xmax=max(stop))) +
            ggtitle(paste0("Window-based \u0394AIC between ", long_wsize, " and ", short_wsize, "bp")) +
            xlab("Chromosomal position (bp)") + ylab("\u0394AIC") +
            geom_rect(aes(fill=bg_color, ymin=-Inf, ymax=Inf, xmin=start, xmax=stop),
                          alpha=0.5, inherit.aes = F) +
            scale_fill_identity() +
            geom_point(shape=16) +
            guides(colour="none", shape="none") +
            theme(
                plot.title = element_text(face = "bold"),
                plot.margin = margin(0.5, 0, 0.5, 0, "cm")
            )
  
  return(plot)
}

# functions: plot multiple trees next to each other
# required packages: ggtree
f_plot_multiple_trees <- function(ls_treefile, ls_annotation, min_branch_support) {
  # set up variables
  plot <- NULL
  df_all <- NULL
  xmax <- 0

  # iterate over treefile
  for (i in 1:length(ls_treefile)) {
    # read tree in Newick format
    tree <- tryCatch(
      {
      treeio::read.newick(ls_treefile[i], node.label='support')
      },
      error = function(x) {
        treeio::read.newick(ls_treefile[i])
      }
    )
    df_tree <- ggtree::fortify(tree)
    
    # set the window size for the tree
    tree_text <- grid::textGrob(ls_annotation[i])
    
    # plot the tree
    if (i == 1) {
      plot <- ggtree(tree, size=1) 
    } else {
      df_tree$x <- df_tree$x + xmax + 0.5
      plot <- plot + geom_tree(data=df_tree, size=1)
    }

    if ("support" %in% colnames(df_tree)) {
      df_tree_subset <- subset(df_tree, support >= min_branch_support)

      if (nrow(df_tree_subset) > 0) {
        plot <- plot + geom_nodelab(aes(label="*"),
                                  data=df_tree_subset,
                                  hjust=1.9, vjust=0.2, size=6)
      }
    }
    
    # add annotation under the tree
    plot <- plot +
      annotation_custom(grob = tree_text,  xmin = min(df_tree$x), xmax = max(df_tree$x), ymin = min(df_tree$y)-0.3, ymax = min(df_tree$y)-0.3)
    
    # add species name for the last tree
    if (i == length(ls_treefile)) {
      plot <- plot + geom_tiplab(data=df_tree, fontface = "italic")
    }
    
    # update variables
    xmax <- max(df_tree$x)
    df_all <- bind_rows(df_all, df_tree)
  }

  # filter trees without labels
  df_all <- subset(df_all, !is.na(label))

  # plot all trees
  plot <- plot +
    geom_line(aes(x, y, color=label), data=df_all, alpha=0.5, linewidth=1) +
    ggtitle("Comparison of window trees") +
    guides(color="none") +
    theme(
      plot.title = element_text(face = "bold"), 
      plot.margin = margin(0.5, 0, 0.5, 0, "cm")
    )

  return(plot)
}

# functions: set colours for each nucleobase
f_extract_base_color <- function(base) {
  switch (base,
    "A" = "red",
    "C" = "blue",
    "G" = "orange",
    "T" = "darkgreen",
    "-" = "lightgrey"
  )
}

# functions: print FASTA alignment as data.frame
# required packages: kableExtra, seqinr
f_print_fasta_alignment <- function(fn_fasta) {
  # read fasta
  fasta_aln <- seqinr::read.fasta(fn_fasta, whole.header = T)

  # convert fasta to data.frame
  df_fasta <- do.call(rbind, fasta_aln)
  df_fasta <- as.data.frame(toupper(df_fasta))
  colnames(df_fasta) <- seq(1, ncol(df_fasta))

  # set up variables
  len_taxa <- nrow(df_fasta)
  len_fasta <- ncol(df_fasta)
  vct_cols <- seq(1, len_fasta)

  # extract polymorphic sites
  idx_polymorphic <- NULL
  for (i in 1:len_fasta) {
    uq_bases <- subset(df_fasta[,i], df_fasta[,i] != "-")
    
    if (length(unique(uq_bases)) > 1) {
      idx_polymorphic <- c(idx_polymorphic, i)
    }
  }

  # extract spaced column names
  df_colnames <- list(" " = 1)
  min_multiplier <- floor(len_fasta/10)
  for (i in 0:min_multiplier) {
    if (i == 0) {
      df_colnames[["1"]] <- 9
    } else if (i == min_multiplier) {
      if (i*10 < len_fasta) {
        df_colnames[[paste(i*10)]] <- len_fasta - (i*10)
      }
      df_colnames[[paste(len_fasta)]] <- 1
    } else {
      df_colnames[[paste(i*10)]] <- 10
    }
  }

  # visualization
  output <- df_fasta %>%
    mutate(across(all_of(vct_cols[!vct_cols %in% idx_polymorphic]), function(z) {
      cell_spec(z, color = "lightgrey")
    })) %>%
    mutate(across(all_of(idx_polymorphic), function(x) {
      sapply(x, function(y) {
        cell_spec(y, bold = T, color = f_extract_base_color(y))
      })
    })) %>%
    kable(escape = F, col.names=NULL) %>%
    kable_styling(full_width = T) %>%
    add_header_above(df_colnames, align="l") %>%
    scroll_box(width = "100%", fixed_thead = TRUE, height = ifelse(len_taxa > 7, "300px", ""))

  return(output)
}