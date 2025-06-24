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