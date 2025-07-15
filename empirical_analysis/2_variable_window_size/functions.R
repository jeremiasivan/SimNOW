# functions for codes_empirical/2_variable_window_size

# function: create data.frame for all windows
f_create_perwindow_df_v2 <- function(numw, fasta_len, window_size) {
  iterator <- 1:numw
  width <- ceiling(log(numw) / log(10)) + 1

  # set up start and end positions
  start_pos <- seq(1, fasta_len - window_size + 1, by = window_size)
  end_pos <- start_pos + window_size - 1

  # add the last window when it is not divisible by the window size
  if (max(end_pos) != fasta_len) {
    start_pos <- c(start_pos, max(end_pos)+1)
    end_pos <- c(end_pos, fasta_len)
  }

  # create a data.frame
  df_windows <- data.frame(window_name = paste0("window_", formatC(iterator, width=width, format="d", flag="0")),
                           start = start_pos, end = end_pos)
  df_windows$length <- df_windows$end - df_windows$start + 1
  
  return(df_windows)
}

# function: plot the topologies across chromosomes
f_plot_top_pos <- function(df, colour_scheme, fn_output) {
  plot <- ggplot(df) +
    geom_rect(aes(xmin=start, xmax=stop, ymin=y-0.6, ymax=y+0.6, fill=topology)) +
    labs(x="Position", y="Chromosome", color="Topology") +
    scale_y_continuous(breaks=unique(df$y), labels=unique(df$chr), expand=c(0.02, 0.02)) +
    theme(axis.title.x=element_text(size=40, margin = margin(t=20, r=0, b=0, l=0)),
          axis.title.y=element_text(size=40, margin = margin(t=0, r=20, b=0, l=0)),
          axis.text.y=element_text(size=40),
          axis.text.x=element_text(size=40),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          legend.title=element_text(size=30),
          legend.text=element_text(size=30),
          legend.key.size=unit(2,"cm"))

  # update the colour scheme
  if (!is.null(colour_scheme)) {
      plot <- plot + scale_fill_manual(values=colour_scheme)
  }
    
  # save the plot
  tiff(filename=fn_output, units="px", width=3840, height=1080)
  print(plot)
  dev.off()
}