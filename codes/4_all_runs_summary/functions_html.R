# functions for codes/4_all_runs_summary

# function: plot accuracy vs window size
# required library: ggplot2
f_acc_wsize_2 <- function(df, type) {
  plot <- NULL
  if (type == "site") {
    plot <- ggplot(df, aes(y=accuracy, x=window_size, ymin=0, ymax=100)) +
      ggtitle("Site Accuracy of Non-Overlapping Windows") +
      xlab("Window size (kb)") + ylab("Accuracy (%)")
  } else if (type == "rmse") {
    plot <- ggplot(df, aes(y=rmse, x=window_size, ymin=0)) +
      ggtitle("RMSE of Non-Overlapping Windows") +
      xlab("Window size (kb)") + ylab("RMSE (%)")
  }

  plot <- plot +
    geom_boxplot(outlier.shape = NA) +
    geom_point(shape=16, aes(colour=r)) +
    facet_wrap(.~r, scales = "free", ncol = 1) +
    guides(colour="none") +
    theme(
      plot.title = element_text(face = "bold"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )

  return(plot)
}

# function: plot site accuracy vs RMSE
# required library: ggplot2
f_acc_rmse_2 <- function(df) {
  plot <- ggplot(df, aes(y=log10(accuracy), x=log10(rmse), ymin=0, xmin=0)) +
    ggtitle("Correlation between Site Accuracy and RMSE") +
    xlab("Log10-transformed RMSE") + ylab("Log10-transformed accuracy") +
    geom_point(aes(colour=r)) +
    geom_smooth(method = lm) +
    facet_wrap(.~r, scales = "free") +
    guides(colour="none") +
    theme(
      plot.title = element_text(face = "bold"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )

  return(plot)
}

# function: plot accuracy vs ic
# required library: ggplot2, viridis
f_ic_acc_2 <- function(df, i, ifc, type) {
  plot <- NULL
  if (type == "site") {
    plot <- ggplot(df, aes(y=accuracy, x=as.numeric(unlist(df[,..i])), ymin=0, ymax=100, group=simulation)) +
      ggtitle(paste("Correlation between", toupper(ifc), "and Site Accuracy")) +
      xlab(toupper(ifc)) + ylab("Accuracy (%)")
  } else if (type == "rmse") {
    plot <- ggplot(df, aes(y=rmse, x=as.numeric(unlist(df[,..i])), ymin=0, group=simulation)) +
      ggtitle(paste("Correlation between", toupper(ifc), "and RMSE")) +
      xlab(toupper(ifc)) + ylab("RMSE (%)")
  }
  
  plot <- plot +
    geom_line(aes(alpha=0.2)) +
    geom_point(aes(colour=window_size)) +
    facet_wrap(.~r, scales = "free") +
    viridis::scale_color_viridis(discrete = TRUE) + 
    labs(color="wsize (kb)") +
    guides(alpha="none") +
    theme(
      plot.title = element_text(face = "bold"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )
  
  return(plot)
}

# function: plot accuracy vs delta ic
# required library: dplyr, ggplot2, viridis
f_delta_ic_acc_2 <- function(df, i, ifc, type) {
  col_ifc <- sym(ifc)
  
  sdf <- df %>%
    group_by(simulation) %>%
    mutate(min_ifc = min(!!col_ifc)) %>%
    mutate(d_ifc = abs(!!col_ifc-min_ifc))
  
  plot <- NULL
  if (type == "site") {
    plot <- ggplot(sdf, aes(y=accuracy, x=d_ifc, ymin=0, ymax=100, group=simulation)) +
      ggtitle(paste0("Correlation between \u0394", toupper(ifc), " and Site Accuracy")) +
      xlab(paste0("\u0394", toupper(ifc))) + ylab("Accuracy (%)")
  } else if (type == "rmse") {
    plot <- ggplot(sdf, aes(y=rmse, x=d_ifc, ymin=0, group=simulation)) +
      ggtitle(paste0("Correlation between \u0394", toupper(ifc), " and RMSE")) +
      xlab(paste0("\u0394", toupper(ifc))) + ylab("RMSE (%)")
  }
  
  plot <- plot +
    geom_line(aes(alpha=0.2)) +
    geom_point(aes(colour=window_size)) +
    facet_wrap(.~r, scales = "free") +
    viridis::scale_color_viridis(discrete = TRUE) + 
    labs(color="wsize (kb)") +
    guides(alpha="none") +
    theme(
      plot.title = element_text(face = "bold"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )
  
  return(plot)
}

# function: plot accuracy loss
# required library: ggplot2
f_acc_loss_2 <- function(df, ifc, type) {
  plot <- NULL
  if (type == "site") {
    plot <- ggplot(df, aes(y=site_loss, x=factor(r), ymin=0)) +
      ggtitle(paste("Loss of Site Accuracy for Lowest", toupper(ifc))) +
      xlab("Recombination rate (rho)") + ylab("Loss of Site Accuracy (%)")
  } else if (type == "rmse") {
    plot <- ggplot(df, aes(y=rmse_gain, x=factor(r), ymin=0)) +
      ggtitle(paste("Increase in RMSE for Lowest", toupper(ifc))) +
      xlab("Recombination rate (rho)") + ylab("Increase in RMSE (%)")
  }
  
  plot <- plot +
    geom_boxplot(outlier.shape = NA) +
    geom_point(shape=16, aes(colour=r)) +
    guides(colour="none") +
    theme(
      plot.title = element_text(face = "bold"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )
  
  return(plot)
}

# function: plot window size vs ic
# required library: ggplot2, viridis
f_ic_wsize_2 <- function(df, i, ifc) {
  plot <- ggplot(df, aes(x=window_size, y=as.numeric(unlist(df[,..i])), group=simulation)) +
    geom_line(aes(alpha=0.2)) +
    geom_point(aes(colour=window_size)) +
    facet_wrap(.~r, scales = "free", ncol = 1) +
    viridis::scale_color_viridis(discrete = TRUE) +
    ggtitle(paste("Correlation between", toupper(ifc), "and Window Size")) + ylab(toupper(ifc)) + xlab("Window Size (kb)") +
    guides(color="none", alpha="none") +
    theme(
      plot.title = element_text(face = "bold"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )
  
  return(plot)
}

# function: plot ranked accuracy vs ranked ic
# required library: ggplot2, viridis
f_rank_ic_acc_2 <- function(df, i, ifc, type) {
  plot <- NULL
  if (type == "site") {
    plot <- ggplot(df, aes(y=accuracy, x=unlist(df[,i]), ymin=0, group=simulation)) +
      ggtitle(paste("Correlation between Ranked", toupper(ifc), "and Ranked Site Accuracy")) +
      xlab(paste("Ranked",toupper(ifc))) + ylab("Ranked Accuracy")
  } else if (type == "rmse") {
    plot <- ggplot(df, aes(y=rmse, x=unlist(df[,i]), ymin=0, group=simulation)) +
      ggtitle(paste("Correlation between Ranked", toupper(ifc), "and Ranked RMSE")) +
      xlab(paste("Ranked",toupper(ifc))) + ylab("Ranked RMSE")
  }
  
  plot <- plot +
    geom_line(aes(alpha=0.2)) +
    geom_point(aes(colour=window_size)) +
    facet_wrap(.~r, scales = "free") +
    viridis::scale_color_viridis(discrete = TRUE) +
    labs(color="wsize (kb)") +
    guides(alpha="none") +
    theme(
      plot.title = element_text(face = "bold"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )
  
  return(plot)
}

# function: plot all plots
f_plot_all <- function(simsum, acc_loss, rank_simsum, i, ifc) {
  # information criteria vs accuracy / rmse
  cat('#### Accuracy {.tabset .tabset-pills} \n')
  cat('##### Site accuracy \n')
  plot(f_ic_acc_2(simsum, i, ifc, "site"))
  cat('  \n\n')
  cat('##### RMSE \n')
  plot(f_ic_acc_2(simsum, i, ifc, "rmse"))
  cat('  \n\n')

  # delta information criteria vs delta accuracy / rmse
  cat('#### Delta accuracy {.tabset .tabset-pills} \n')
  cat('##### Site accuracy \n')
  plot(f_delta_ic_acc_2(simsum, i, ifc, "site"))
  cat('  \n\n')
  cat('##### RMSE \n')
  plot(f_delta_ic_acc_2(simsum, i, ifc, "rmse"))
  cat('  \n\n')
  
  # accuracy loss
  cat('#### Accuracy loss {.tabset .tabset-pills} \n')
  cat('##### Site accuracy \n')
  plot(f_acc_loss_2(acc_loss, ifc, "site"))
  cat('  \n\n')
  cat('##### RMSE \n')
  plot(f_acc_loss_2(acc_loss, ifc, "rmse"))
  cat('  \n\n')
  
  # ranked IC vs ranked accuracy
  cat('#### Ranked accuracy {.tabset .tabset-pills} \n')
  cat('##### Site accuracy \n')
  plot(f_rank_ic_acc_2(rank_simsum, i, ifc, "site"))
  cat('  \n\n')
  cat('##### RMSE \n')
  plot(f_rank_ic_acc_2(rank_simsum, i, ifc, "rmse"))
  cat('  \n\n')

  # information criterion vs window size
  cat('#### Window size \n')
  plot(f_ic_wsize_2(simsum, i, ifc))
  cat('  \n\n')
}