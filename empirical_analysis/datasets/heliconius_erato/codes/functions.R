# function: calculate topology distribution
f_count_topology_freq <- function(window_size, list_chr, step, outdir, fn_output) {
  df_output <- data.table::data.table(topology=character(), n=numeric(), cum.percentage=numeric())

  for (c in list_chr) {
    df_uqtops <- data.table::fread(paste0(outdir, "/", c, "/", step, "/", window_size, "/", window_size, ".uqtops"))
    df_output <- rbind(df_output, df_uqtops)
  }

  df_output <- df_output %>%
    select(topology, n) %>%
    group_by(topology) %>%
    summarise_all(sum) %>%
    arrange(desc(n)) %>%
    mutate(cum.percentage=cumsum(n)/sum(n)*100)

  data.table::fwrite(df_output, file=fn_output, quote=F, sep="\t")
}

# function: calculate topology weights
f_compare_topology_weight <- function(df_input, df_uqtops, colname) {
  if (nrow(df_input) == 0) {
    # remove the last column
    df_uqtops$cum.percentage <- NULL

    # update colnames
    colnames(df_uqtops) <- c("topology", colname)
    return(df_uqtops)
  }

  # iterate over unique topology
  top_dist <- c()
  for (r in 1:nrow(df_input)) {
    idx <- match(df_input$topology[r], df_uqtops$topology)
    if (is.na(idx)) {
      top_dist <- c(top_dist, 0)
      next
    }

    top_dist <- c(top_dist, df_uqtops$n[idx])
  }

  # add new topologies
  new_tops <- subset(df_uqtops$topology, !df_uqtops$topology %in% df_input$topology)
  for (r in new_tops) {
    # new entry
    new_row <- as.data.frame(t(c(r, rep(0, ncol(df_input)-1))))
    df_input <- rbind(df_input, new_row, use.names=FALSE)

    idx <- match(r, df_uqtops$topology)
    top_dist <- c(top_dist, df_uqtops$n[idx])
  }

  # combine the two data.table
  df_input[, (colname) := top_dist]
  return(df_input)
}