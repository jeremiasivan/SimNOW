#################################
#  Two-window HiMAST from NOW   #
#################################

fn_mstops <- "~/sim.mstops"
fn_atsum <- "~/sim.atsum"

dir_fasta <- "~/sim/windows/alignment"
dir_iqtree2 <- "~/iqtree-2.2.5.hmmster-Linux/bin/iqtree2"
dir_seqkit <- "~/seqkit"

dir_output <- "~/simulation"

#################################
# install.packages("data.table")
# install.packages("log4r")
library(dplyr)

# logger
fn_log <- paste0(dir_output, "/log.txt")
log_appender <- log4r::file_appender(fn_log, append = TRUE, layout = log4r::default_log_layout())
fn_logger <- log4r::logger(threshold = "INFO", appenders = log_appender)

# output directory for HiMAST-based windows
dir_output_windows <- paste0(dir_output, "/windows/")
if (!dir.exists(dir_output_windows)) {
  dir.create(dir_output_windows, recursive=T)
}

# output directory for concatenated sequences
dir_output_runs <- paste0(dir_output, "/runs/")
if (!dir.exists(dir_output_runs)) {
  dir.create(dir_output_runs, recursive=T)
}

# temporary output directory
dir_output_tmp <- paste0(dir_output, "/tmp/")
if (!dir.exists(dir_output_tmp)) {
  dir.create(dir_output_tmp, recursive=T)
}

#################################

# read atsum file
atsum <- data.table::fread(fn_atsum)

# output data.table
dt_output <- data.table::data.table(source = character(),
                                    start = numeric(),
                                    stop = numeric(),
                                    length = numeric(),
                                    topology = character())

# initialize variables
current_topology <- ""
current_length <- ""
current_start <- 0

p <- nrow(atsum)
options(scipen=9999)

#################################

# functions
run_seqkit <- function(seq1, seq2, seqkit, output) {
  system(paste0(seqkit, " concat ", seq1, " ", seq2, " > ", output))
}

run_himast <- function(seq, top, iqtree2) {
  mast_cl <- paste(iqtree2,
                    "-s", seq,
                    "-m 'JC+T'",
                    "-te", top,
                    "-hmmster -hmm_min_stran 0.99999 -T 1 -redo")
  system(mast_cl)
}

extract_hmmtrees <- function(fn_hmm, fn_logger, fn_aln_in, fn_atsum1, fn_atsum2, dt_output, ls_topology, current_start, dir_win, dir_tmp) {
  # extract HMM per-site topology
  hmmtrees <- system(paste("grep '^\\['", fn_hmm), intern=TRUE)
  trees <- data.table::data.table(do.call('rbind', strsplit(as.character(hmmtrees), '\\[|,|\\]')))
  trees[,1] <- NULL
  
  colnames(trees) <- c("start","stop","tree")
  trees$start <- as.numeric(trees$start)
  trees$stop <- as.numeric(trees$stop)
  trees$tree <- as.numeric(trees$tree)
  
  # check if the first block is the first tree
  if (trees$tree[1] != 1) {
    log4r::warn(fn_logger, paste("Warn: first block for step", i, "is not first topology."))
  }
  
  # check the number of hmm partitions
  n <- nrow(trees)
  if (n > 2) {
    log4r::warn(fn_logger, paste0("Warn: ", n, " partitions for step ", i, '.'))
  }

  # update variables
  current_length <- 0
  current_topology <- ""

  if (n == 1) {
    # copy the file to other folder
    system(paste0("cp ", fn_aln_in, " ", dir_tmp, fn_atsum2))
    
    # update variables
    current_length <- trees$stop[1] - trees$start[1] + 1
    current_topology <- ls_topology[trees$tree[1]]
  
  } else {
    # read input FASTA file
    s <- seqinr::read.fasta(fn_aln_in, whole.header = T)
    
    # loop over hmm partitions
    for (j in 1:n) {
      # extract the block
      subfasta <- lapply(s, function(x) x[seq(from = trees$start[j], to = trees$stop[j])])
      subfasta_df <- do.call(rbind,subfasta)
      subfasta <- setNames(split(subfasta_df, seq(nrow(subfasta_df))), rownames(subfasta_df))
      
      # save the sequence in other file
      fn_outfile <- ""
      if (j != n) {
        fn_fasta <- strsplit(fn_atsum1, ".fa")[[1]]
        fn_outfile <- paste0(dir_win, fn_fasta, "_", j, ".fa")
        tmp_length <- trees$stop[j] - trees$start[j] + 1
        
        # add entry to output data.table
        dt_output <- rbind(dt_output, list(paste0(fn_fasta, "_", j, ".fa"),
                                        current_start+1,
                                        current_start+tmp_length,
                                        tmp_length,
                                        ls_topology[trees$tree[j]]))
        
        # update variables
        current_start <- current_start + tmp_length
        
      } else {
        fn_outfile <- paste0(dir_tmp, fn_atsum2)
        
        # update variables
        current_length <- trees$stop[j] - trees$start[j] + 1
        current_topology <- ls_topology[trees$tree[j]]
      }
      
      # write the sequence
      seqinr::write.fasta(sequences=subfasta, names=names(subfasta), file.out=fn_outfile, nbchar=100)
    }
  }

  return(list(out=dt_output, length=current_length, topology=current_topology, start=current_start))
}

#################################

# run
for (i in 1:p) {
  fn_current_fasta <- paste0(dir_output_tmp, atsum$source[i])
  
  # check the window from previous run 
  if (file.exists(fn_current_fasta)) {
    # save the last window and check if the windows have the same topology
    if (i == p || current_topology == atsum$topology[i+1]) {
      # copy the file to other folder
      system(paste0("cp ", fn_current_fasta, " ", dir_output_windows))
      
      # remove the file from temporary folder
      unlink(fn_current_fasta)
      
      # add entry to output data.table
      dt_output <- rbind(dt_output, list(atsum$source[i],
                                      current_start+1,
                                      current_start+current_length,
                                      current_length,
                                      current_topology))
      
      # update variables
      current_start <- current_start + current_length
      next
    
    } else {
      # create output directory for the run
      dir_run <- paste0(dir_output_runs, i, "/")
      if (!dir.exists(dir_run)) {
        dir.create(dir_run)
      }
      
      # run seqkit to concat the sequences
      seq1 <- paste0(dir_output_tmp, "/", atsum$source[i])
      seq2 <- paste0(dir_fasta, "/", atsum$source[i+1])
      fn_alignment <- paste0(dir_run, "concat_", i, ".fa")
      run_seqkit(seq1, seq2, dir_seqkit, fn_alignment)
      
      # topology file
      fn_topology <- paste0(dir_run, paste0("topology_", i, ".txt"))
      ls_topology <- c(current_topology, atsum$topology[i+1])
      write.table(ls_topology, file = fn_topology, quote = F, row.names = F, col.names = F)
      
      # run HiMAST
      run_himast(fn_alignment, fn_topology, dir_iqtree2)
      
      # check if HMM file exists
      fn_hmm <- paste0(fn_alignment, ".hmm")
      if (!file.exists(fn_hmm)) {
        errmsg <- paste("Error: HMM file for step", i, "is not found. Exited.")

        log4r::error(fn_logger, errmsg)
        stop(errmsg)
      }

      # extract HMM topologies
      output <- extract_hmmtrees(fn_hmm, fn_logger,
                                fn_alignment, atsum$source[i], atsum$source[i+1],
                                dt_output, ls_topology,
                                current_start,
                                dir_output_windows, dir_output_tmp)

      dt_output <- output$out
      current_length <- output$length
      current_topology <- output$topology
      current_start <- output$start
    }
    
  } else {
    # save the last window and check if the windows have the same topology
    if (i == p || atsum$topology[i] == atsum$topology[i+1]) {
      # copy the file to other folder
      system(paste0("cp ", dir_fasta, "/", atsum$source[i], " ", dir_output_windows))
      
      # add entry to output data.table
      dt_output <- rbind(dt_output, atsum[i,])
      
      # update variables
      current_start <- current_start + atsum$length[i]
      next
    
    } else {
      # create output directory for the run
      dir_run <- paste0(dir_output_runs, i, "/")
      if (!dir.exists(dir_run)) {
        dir.create(dir_run)
      }
      
      # run seqkit to concat the sequences
      seq1 <- paste0(dir_fasta, "/", atsum$source[i])
      seq2 <- paste0(dir_fasta, "/", atsum$source[i+1])
      fn_alignment <- paste0(dir_run, "concat_", i, ".fa")
      run_seqkit(seq1, seq2, dir_seqkit, fn_alignment)
      
      # topology file
      fn_topology <- paste0(dir_run, paste0("topology_", i, ".txt"))
      ls_topology <- c(atsum$topology[i], atsum$topology[i+1])
      write.table(ls_topology, file = fn_topology, quote = F, row.names = F, col.names = F)
      
      # run HiMAST
      run_himast(fn_alignment, fn_topology, dir_iqtree2)
      
      # check if HMM file exists
      fn_hmm <- paste0(fn_alignment, ".hmm")
      if (!file.exists(fn_hmm)) {
        errmsg <- paste("Error: HMM file for step", i, "is not found. Exited.")

        log4r::error(fn_logger, errmsg)
        stop(errmsg)
      }
      
      # extract HMM per-site topology
      output <- extract_hmmtrees(fn_hmm, fn_logger,
                                fn_alignment, atsum$source[i], atsum$source[i+1],
                                dt_output, ls_topology,
                                current_start,
                                dir_output_windows, dir_output_tmp)

      dt_output <- output$out
      current_length <- output$length
      current_topology <- output$topology
      current_start <- output$start
    }
  }
}

# write the output file
data.table::fwrite(dt_output, paste0(dir_output, "/output.atsum"), quote=F, sep="\t")

#################################

# summarize the accuracy
msfile <- data.table::fread(fn_mstops)

mstree_uq <- unique(msfile$topology)
mstree_uq_ln <- length(mstree_uq)
smtree_uq <- unique(dt_output$topology)

# match window topology with ms topology
for (j in mstree_uq) {
  for (k in smtree_uq) {
    if (ape::all.equal.phylo(ape::unroot(ape::read.tree(text = j)), ape::read.tree(text = k))) {
      dt_output[, topology := gsub(k, j, topology, fixed = TRUE)]
      break
    }
  }
}

# create empty dataframe
output <- data.table::data.table(c(mstree_uq,"NT"), matrix(0, nrow=mstree_uq_ln+1, ncol=mstree_uq_ln))
data.table::setnames(output, c("topology", mstree_uq))

# iterate through ms topology file
for (k in 1:nrow(msfile)){
  # extract the hmm topology based on ms topology switching
  start_idx <- tail(which(dt_output$start <= msfile$start[k]), 1)
  stop_idx <- head(which(dt_output$stop >= msfile$stop[k]), 1)
  
  # update lengths
  sites <- dt_output[start_idx:stop_idx,]
  sites$start[1] <- msfile$start[k]
  sites$length[1] <- sites$stop[1] - sites$start[1] + 1
  
  last_idx <- nrow(sites)
  sites$stop[last_idx] <- msfile$stop[k]
  sites$length[last_idx] <- sites$stop[last_idx] - sites$start[last_idx] + 1
  
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
  
  rm(sites)
}

data.table::fwrite(output, file=paste0(dir_output, "/output.cmp"), quote=F, sep="\t")

# calculate accuracy
output <- as.matrix(output[,-1])
acc <- sum(as.numeric(diag(output))) / sum(as.numeric(output)) * 100
acc <- round(acc, 3)

log4r::info(fn_logger, paste("Accuracy:",acc))
