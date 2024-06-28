# function: plot the topology distribution for NOW
f_topsimsum <- function(df) {
  row_sums <- sum(df$count)

  df <- df %>%
    group_by(topology, colour) %>%
    summarise(count=sum(count)) %>%
    mutate(cum.percentage = round(count/row_sums*100, 3)) %>%
    arrange(desc(count))

  df$topology <- forcats::fct_inorder(df$topology)

  plot1 <- ggplot(df, aes(x=topology, y=count))
  plot2 <- NULL
  
  # check if colour column is empty
  if (!all(is.na(df$colour))) {
    plot1 <- plot1 + geom_bar(position="dodge", stat="identity", aes(fill=colour)) + scale_fill_identity()
    plot2 <- ggplot(subset(df, colour!="#000000"), aes(x=topology, y=count)) + geom_bar(position="dodge", stat="identity", aes(fill=colour)) + scale_fill_identity()
  } else {
    plot1 <- plot1 + geom_bar(position="dodge", stat="identity", aes(fill=topology)) + viridis::scale_fill_viridis(discrete = TRUE)
    plot2 <- ggplot(df[1:8,], aes(x=topology, y=count)) + geom_bar(position="dodge", stat="identity", aes(fill=topology)) + viridis::scale_fill_viridis(discrete = TRUE)
  }

  plot1 <- plot1 + 
    geom_text(aes(label=cum.percentage), position=position_dodge(width=0.9), vjust=-0.5, size=5) +
    ggtitle("Topology distribution of erato-sara Heliconius genome") + ylab("Count") + xlab("Topology") +
    guides(fill="none") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 50),
      plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),
      axis.title.x=element_text(size = 40, margin = margin(t=20, r=0, b=0, l=0)),
      axis.title.y=element_text(size = 40, margin = margin(t=0, r=20, b=0, l=0)),
      axis.text.y=element_text(size=30),
      axis.text.x=element_text(size=30, angle = 90, vjust = 0.5, hjust=1),
      legend.title=element_text(size=30),
      legend.text=element_text(size=30),
      legend.key.size=unit(2,"cm")
    )

  plot2 <- plot2 +
    geom_text(aes(label=cum.percentage), position=position_dodge(width=0.9), vjust=-0.5, size=10) +
    ggtitle("Top 8 topologies of erato-sara Heliconius genome") + ylab("Count") + xlab("Topology") +
    guides(fill="none") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 50),
      plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),
      axis.title.x=element_text(size = 40, margin = margin(t=20, r=0, b=0, l=0)),
      axis.title.y=element_text(size = 40, margin = margin(t=0, r=20, b=0, l=0)),
      axis.text.y=element_text(size=30),
      axis.text.x=element_text(size=30, angle = 90, vjust = 0.5, hjust=1),
      legend.title=element_text(size=30),
      legend.text=element_text(size=30),
      legend.key.size=unit(2,"cm")
    )

  return(list(plot1=plot1, plot2=plot2, df=df))
}

# function: plot the consecutive topology distribution for NOW
f_cnsimsum <- function(df1, df2) {
  df1 <- df1 %>%
    group_by(topology, count, colour, tree) %>%
    summarise(n=n())

  df1$tree <- factor(df1$tree, levels=stringr::str_sort(unique(df1$tree), numeric=TRUE))

  plot1 <- ggplot(df1, aes(x=count, y=n))
  plot2 <- NULL
  
  # check if colour column is empty
  if (!all(is.na(df1$colour))) {
    plot1 <- plot1 + geom_bar(position="dodge", stat="identity", aes(fill=colour)) + scale_fill_identity()
    plot2 <- ggplot(subset(df1, colour!="#000000"), aes(x=count, y=n)) + geom_bar(position="dodge", stat="identity", aes(fill=colour)) + scale_fill_identity()
  } else {
    top_eight_tops <- df2$topology[1:8]
    df2 <- subset(df1, topology %in% top_eight_tops)

    plot1 <- plot1 + geom_bar(position="dodge", stat="identity", aes(fill=topology)) + viridis::scale_fill_viridis(discrete = TRUE)
    plot2 <- ggplot(df2, aes(x=count, y=n)) + geom_bar(position="dodge", stat="identity", aes(fill=topology)) + viridis::scale_fill_viridis(discrete = TRUE)
  }

  plot1 <- plot1 +
    facet_wrap(.~topology, scales = "free") +
    ggtitle("Number of consecutive windows per topology") + ylab("Frequency") + xlab("Number of consecutive windows") +
    guides(fill="none") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 50),
      plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),
      axis.title.x=element_text(size = 40, margin = margin(t=20, r=0, b=0, l=0)),
      axis.title.y=element_text(size = 40, margin = margin(t=0, r=20, b=0, l=0)),
      axis.text.y=element_text(size=30),
      axis.text.x=element_text(size=30),
      legend.title=element_text(size=30),
      legend.text=element_text(size=30),
      legend.key.size=unit(2,"cm")
    )

  plot2 <- plot2 + 
    facet_wrap(.~tree, scales = "free") +
    ggtitle("Number of consecutive windows for 8 top topologies") + ylab("Frequency") + xlab("Number of consecutive windows") +
    guides(fill="none") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 50),
      plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),
      axis.title.x=element_text(size = 40, margin = margin(t=20, r=0, b=0, l=0)),
      axis.title.y=element_text(size = 40, margin = margin(t=0, r=20, b=0, l=0)),
      axis.text.y=element_text(size=30),
      axis.text.x=element_text(size=30),
      strip.text = element_text(size = 20),
      legend.title=element_text(size=30),
      legend.text=element_text(size=30),
      legend.key.size=unit(2,"cm")
  )

  return(list(plot1=plot1, plot2=plot2))
}

# function: calculate topology distribution for stepwise NOW
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

# function: calculate topology weights for stepwise NOW
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