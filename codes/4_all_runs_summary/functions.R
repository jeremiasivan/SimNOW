# functions for codes/4_all_runs_summary

# function: plot accuracy vs window size
# required library: ggplot2
f_acc_wsize <- function(df, type) {
  plot <- NULL
  if (type == "site") {
    plot <- ggplot(df, aes(y=accuracy, x=window_size, ymin=0, ymax=100)) +
      ggtitle("Site Classification Accuracy of Non-Overlapping Windows in Simulated Datasets") +
      xlab("Window size (kb)") + ylab("Accuracy (%)")
  } else if (type == "rmse") {
    plot <- ggplot(df, aes(y=rmse, x=window_size, ymin=0)) +
      ggtitle("RMSE of Non-Overlapping Windows in Simulated Datasets") +
      xlab("Window size (kb)") + ylab("RMSE (%)")
  }

  plot <- plot +
    geom_boxplot(outlier.shape = NA) +
    geom_point(shape=16, aes(colour=r, size=20)) +
    facet_wrap(.~r, scales = "free", ncol = 1) +
    guides(size="none", colour="none") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 50),
      plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),
      axis.title.x=element_text(size = 40, margin = margin(t=20, r=0, b=0, l=0)),
      axis.title.y=element_text(size = 40, margin = margin(t=0, r=20, b=0, l=0)),
      axis.text.y=element_text(size=30),
      axis.text.x=element_text(size=30),
      strip.text = element_text(size = 30)
    )

  return(plot)
}

# function: plot site accuracy vs RMSE
# required library: ggplot2
f_acc_rmse <- function(df) {
  plot <- ggplot(df, aes(y=log10(accuracy), x=log10(rmse), ymin=0, xmin=0)) +
    ggtitle("Correlation between Site Classification Accuracy and RMSE") +
    xlab("Log10-transformed RMSE") + ylab("Log10-transformed accuracy") +
    geom_point(aes(colour=r, size=20)) +
    geom_smooth(method = lm) +
    facet_wrap(.~r, scales = "free") +
    guides(size="none", colour="none") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 50),
      plot.margin = margin(1.25, 0, 1.25, 0, "cm"),
      axis.title.x=element_text(size = 40, margin = margin(t=20, r=0, b=0, l=0)),
      axis.title.y=element_text(size = 40, margin = margin(t=0, r=20, b=0, l=0)),
      axis.text.y=element_text(size=30),
      axis.text.x=element_text(size=30),
      strip.text = element_text(size = 30),
      legend.title=element_text(size=30),
      legend.text=element_text(size=30),
      legend.key.size=unit(2,"cm")
    )

  return(plot)
}

# function: plot accuracy vs ic
# required library: ggplot2, viridis
f_ic_acc <- function(df, i, ifc, type) {
  plot <- NULL
  if (type == "site") {
    plot <- ggplot(df, aes(y=accuracy, x=as.numeric(unlist(df[,..i])), ymin=0, ymax=100, group=simulation)) +
      ggtitle(paste("Correlation between", toupper(ifc), "and Site Classification Accuracy")) +
      xlab(toupper(ifc)) + ylab("Accuracy (%)")
  } else if (type == "rmse") {
    plot <- ggplot(df, aes(y=rmse, x=as.numeric(unlist(df[,..i])), ymin=0, group=simulation)) +
      ggtitle(paste("Correlation between", toupper(ifc), "and RMSE")) +
      xlab(toupper(ifc)) + ylab("RMSE (%)")
  }
  
  plot <- plot +
    geom_line(aes(size=1, alpha=0.2)) +
    geom_point(aes(colour=window_size, size=20)) +
    facet_wrap(.~r, scales = "free") +
    viridis::scale_color_viridis(discrete = TRUE) + 
    labs(color="wsize (kb)") +
    guides(alpha="none", size="none", colour=guide_legend(override.aes=list(size=8))) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 50),
      plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),
      axis.title.x=element_text(size = 40, margin = margin(t=20, r=0, b=0, l=0)),
      axis.title.y=element_text(size = 40, margin = margin(t=0, r=20, b=0, l=0)),
      axis.text.y=element_text(size=30),
      axis.text.x=element_text(size=30),
      strip.text = element_text(size = 30),
      legend.title=element_text(size=30),
      legend.text=element_text(size=30),
      legend.key.size=unit(2,"cm")
    )
  
  return(plot)
}

# function: plot accuracy vs delta ic
# required library: dplyr, ggplot2, viridis
f_delta_ic_acc <- function(df, i, ifc, type) {
  col_ifc <- sym(ifc)
  
  sdf <- df %>%
    group_by(simulation) %>%
    mutate(min_ifc = min(!!col_ifc)) %>%
    mutate(d_ifc = abs(!!col_ifc-min_ifc))
  
  plot <- NULL
  if (type == "site") {
    plot <- ggplot(sdf, aes(y=accuracy, x=d_ifc, ymin=0, ymax=100, group=simulation)) +
      ggtitle(paste0("Correlation between \u0394", toupper(ifc), " and Site Classification Accuracy")) +
      xlab(paste0("\u0394", toupper(ifc))) + ylab("Accuracy (%)")
  } else if (type == "rmse") {
    plot <- ggplot(sdf, aes(y=rmse, x=d_ifc, ymin=0, group=simulation)) +
      ggtitle(paste0("Correlation between \u0394", toupper(ifc), " and RMSE")) +
      xlab(paste0("\u0394", toupper(ifc))) + ylab("RMSE (%)")
  }
  
  plot <- plot +
    geom_line(aes(size=1, alpha=0.2)) +
    geom_point(aes(colour=window_size, size=20)) +
    facet_wrap(.~r, scales = "free") +
    viridis::scale_color_viridis(discrete = TRUE) + 
    labs(color="wsize (kb)") +
    guides(alpha="none", size="none", colour=guide_legend(override.aes=list(size=8))) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 50),
      plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),
      axis.title.x=element_text(size = 40, margin = margin(t=20, r=0, b=0, l=0)),
      axis.title.y=element_text(size = 40, margin = margin(t=0, r=20, b=0, l=0)),
      axis.text.y=element_text(size=30),
      axis.text.x=element_text(size=30),
      strip.text = element_text(size = 30),
      legend.title=element_text(size=30),
      legend.text=element_text(size=30),
      legend.key.size=unit(2,"cm")
    )
  
  return(plot)
}

# function: plot accuracy loss
# required library: ggplot2
f_acc_loss <- function(df, ifc, type) {
  plot <- NULL
  if (type == "site") {
    plot <- ggplot(df, aes(y=site_loss, x=factor(r), ymin=0)) +
      ggtitle(paste("Loss of Site Classification Accuracy for Lowest", toupper(ifc))) +
      xlab("Recombination rate (rho)") + ylab("Loss of Site Classification Accuracy (%)")
  } else if (type == "rmse") {
    plot <- ggplot(df, aes(y=rmse_gain, x=factor(r), ymin=0)) +
      ggtitle(paste("Increase in RMSE for Lowest", toupper(ifc))) +
      xlab("Recombination rate (rho)") + ylab("Increase in RMSE (%)")
  }
  
  plot <- plot +
    geom_boxplot(outlier.shape = NA) +
    geom_point(shape=16, aes(colour=r, size=20)) +
    guides(size="none", colour="none") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 50),
      plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),
      axis.title.x=element_text(size = 40, margin = margin(t=20, r=0, b=0, l=0)),
      axis.title.y=element_text(size = 40, margin = margin(t=0, r=20, b=0, l=0)),
      axis.text.y=element_text(size=30),
      axis.text.x=element_text(size=30)
    )
  
  return(plot)
}

# function: plot window size vs ic
# required library: ggplot2, viridis
f_ic_wsize <- function(df, i, ifc) {
  plot <- ggplot(df, aes(x=window_size, y=as.numeric(unlist(df[,..i])), group=simulation)) +
    geom_line(aes(size=1, alpha=0.2)) +
    geom_point(aes(colour=window_size, size=20)) +
    facet_wrap(.~r, scales = "free", ncol = 1) +
    viridis::scale_color_viridis(discrete = TRUE) +
    ggtitle(paste("Correlation between", toupper(ifc), "and Window Size")) + ylab(toupper(ifc)) + xlab("Window Size (kb)") +
    guides(color="none", alpha="none", size="none") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 50),
      plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),
      axis.title.x=element_text(size = 40, margin = margin(t=20, r=0, b=0, l=0)),
      axis.title.y=element_text(size = 40, margin = margin(t=0, r=20, b=0, l=0)),
      axis.text.y=element_text(size=30),
      axis.text.x=element_text(size=30),
      strip.text = element_text(size = 30),
      legend.title=element_text(size=30),
      legend.text=element_text(size=30),
      legend.key.size=unit(2,"cm")
    )
  
  return(plot)
}

# function: plot ranked accuracy vs ranked ic
# required library: ggplot2, viridis
f_rank_ic_acc <- function(df, i, ifc, type) {
  plot <- NULL
  if (type == "site") {
    plot <- ggplot(df, aes(y=accuracy, x=unlist(df[,i]), ymin=0, group=simulation)) +
      ggtitle(paste("Correlation between Ranked", toupper(ifc), "and Ranked Site Classification Accuracy")) +
      xlab(paste("Ranked",toupper(ifc))) + ylab("Ranked Accuracy")
  } else if (type == "rmse") {
    plot <- ggplot(df, aes(y=rmse, x=unlist(df[,i]), ymin=0, group=simulation)) +
      ggtitle(paste("Correlation between Ranked", toupper(ifc), "and Ranked RMSE")) +
      xlab(paste("Ranked",toupper(ifc))) + ylab("Ranked RMSE")
  }
  
  plot <- plot +
    geom_line(aes(size=1, alpha=0.2)) +
    geom_point(aes(colour=window_size, size=20)) +
    facet_wrap(.~r, scales = "free") +
    viridis::scale_color_viridis(discrete = TRUE) +
    labs(color="wsize (kb)") +
    guides(alpha="none", size="none", colour=guide_legend(override.aes=list(size=8))) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 50),
      plot.margin = margin(1.25, 0, 1.25, 0, "cm"),
      axis.title.x=element_text(size = 40, margin = margin(t=20, r=0, b=0, l=0)),
      axis.title.y=element_text(size = 40, margin = margin(t=0, r=20, b=0, l=0)),
      axis.text.y=element_text(size=30),
      axis.text.x=element_text(size=30),
      strip.text = element_text(size = 30),
      legend.title=element_text(size=30),
      legend.text=element_text(size=30),
      legend.key.size=unit(2,"cm")
    )
  
  return(plot)
}