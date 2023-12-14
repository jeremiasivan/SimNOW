# functions for codes_empirical/non_overlapping_window

# function: get factors of an integer
f_get_factors <- function(x, min) {
  x <- as.integer(x)
  div <- seq_len(abs(x))
  factors <- div[x %% div == 0L]
  
  return(factors[factors >= min])
}