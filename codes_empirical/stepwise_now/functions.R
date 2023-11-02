# functions for codes_empirical/steowise_now

# function: create window alignment
# required packages: seqinr
f_window_aln <- function(fasta, start, end, dir_perwindow, i, width) {
  subfasta <- lapply(fasta, function(x) x[seq(from = start, to = end)])
  subfasta <- do.call(rbind,subfasta)
  subfasta <- setNames(split(subfasta, seq(nrow(subfasta))), rownames(subfasta))
  
  window_name <- paste0("window_", formatC(i, width=width, format="d", flag="0"))
  dir_window <- paste0(dir_perwindow, "/", window_name, "/")
  if (!dir.exists(dir_window)) {
    dir.create(dir_window, recursive = T)
  }

  fn_out <- paste0(dir_window, window_name, ".fa")
  seqinr::write.fasta(sequences=subfasta, names=names(subfasta), file.out=fn_out, nbchar=100)
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

# function: create per-window aligment and tree
f_perwindow_run <- function(dir_perwindow, wsize, len_window, fasta, dir_iqtree2) {
    # remove all alignments in the folder
    unlink(paste0(dir_perwindow,"*.fa"))
    
    # generate window alignments
    start <- 1
    wi <- ceiling(log(len_window) / log(10)) + 1
    
    for (j in 1:len_window) {
        # create window alignment
        f_window_aln(fasta, start, j*wsize, dir_perwindow, j, wi)

        # update the start position
        start <- start + wsize
    }
}

# function: create per-window aligment and tree
# required packages: ape, data.table, stringr
f_perwindow_sum <- function(dir_perwindow, wsize, len_taxa, min_informative_sites, fn_output) {
    df_output <- data.table::data.table()

    # list all windows
    start <- 1
    dirs <- stringr::str_sort(list.dirs(dir_perwindow, recursive=FALSE, full.names=FALSE), numeric=TRUE)
    
    # iterate over windows
    for (j in 1:length(dirs)) {
        is_informative <- FALSE

        # alignment file
        fn_fasta <- paste0(dir_perwindow, "/", dirs[j], "/", dirs[j], ".fa")

        # output files from IQ-Tree2
        fn_iqtree <- paste0(fn_fasta,".iqtree")
        fn_treefile <- paste0(fn_fasta,".treefile")

        if (file.exists(fn_treefile)) {
            # check the number of indormative sites
            len_informative_sites <- as.numeric(gsub("^.* ", "", system(paste("grep '^Number of parsimony informative sites'",fn_iqtree), intern = T)))

            # check the number of leaf nodes
            tree <- ape::read.tree(file=fn_treefile)
            len_leaf <- length(tree$tip.label)

            if (len_informative_sites >= min_informative_sites && len_leaf == len_taxa) {
                is_informative <- TRUE
            }
        }

        df_output <- rbind(df_output, list(name=dirs[j],start=start,end=j*wsize,is_informative=is_informative))
        start <- start + wsize
    }

    data.table::fwrite(df_output, fn_output, quote=FALSE, sep="\t")
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
            geom_line() +
            geom_point(shape=16) +
            guides(colour="none", shape="none") +
            theme(
                plot.title = element_text(face = "bold"),
                plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
            )
  
  return(plot)
}

# functions: plot multiple trees next to each other
# required packages: ggtree
f_plot_multiple_trees <- function(ls_treefile, ls_annotation, min_bootstrap) {
  # set up variables
  plot <- NULL
  df_all <- NULL
  xmax <- 0

  # iterate over treefile
  for (i in 1:length(ls_treefile)) {
    # read tree in Newick format
    tree <- treeio::read.newick(ls_treefile[i], node.label='support')
    df_tree <- ggtree::fortify(tree)
    
    # set the window size for the tree
    tree_text <- grid::textGrob(ls_annotation[i])
    
    # plot the tree
    if (i == 1) {
      plot <- ggtree(tree, size=1) +
        geom_nodelab(aes(label="*", subset=support>min_bootstrap), hjust=1.9, vjust=0.2, size=6)
    } else {
      df_tree$x <- df_tree$x + xmax + 0.5
      
      plot <- plot + geom_tree(data=df_tree, size=1) +
        geom_nodelab(data=df_tree, aes(label="*", subset=support>min_bootstrap), hjust=1.9, vjust=0.2, size=6)
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
  df_all <- df_all %>% filter(!is.na(label))

  # plot all trees
  plot <- plot +
    geom_line(data=df_all, aes(x, y, color=label), alpha=0.5, linewidth=1) +
    ggplot2::ggtitle("Comparison of Window Trees") +
    guides(color="none") +
    theme(
      plot.title = element_text(face = "bold"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
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
      df_colnames[[paste(i*10)]] <- 10
      df_colnames[[paste(len_fasta)]] <- len_fasta - (i*10) + 1
    } else {
      df_colnames[[paste(i*10)]] <- 10
    }
  }

  # visualization
  plot <- df_fasta %>%
    mutate(across(vct_cols[!vct_cols %in% idx_polymorphic], function(z) {
      cell_spec(z, color = "lightgrey")
    })) %>%
    mutate(across(idx_polymorphic, function(x) {
      sapply(x, function(y) {
        cell_spec(y, bold = T, color = f_extract_base_color(y))
      })
    })) %>%
    kable(escape = F, col.names=NULL) %>%
    kable_styling(full_width = T) %>%
    add_header_above(df_colnames, align="l") %>%
    scroll_box(width = "100%", fixed_thead = TRUE, height = ifelse(len_taxa > 7, "300px", ""))
  
  return(plot)
}