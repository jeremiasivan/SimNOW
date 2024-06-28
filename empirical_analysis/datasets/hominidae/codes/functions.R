# function: extract window topology
f_extract_perwindow_topology <- function(fn_perwindowsum, dir_output_prefix, chr, step, window_size, thread) {
    # open file
    df_perwindowsum <- data.table::fread(fn_perwindowsum)

    # create doSNOW cluster
    nwcl <- makeCluster(thread)
    doSNOW::registerDoSNOW(nwcl)

    # iterate over windows
    df_output_chr <- foreach (window=df_perwindowsum$name, .combine='rbind') %dopar% {
        # check if treefile exists
        fn_tree <- paste0(dir_output_prefix, chr, "/", step, "/", window_size, "/filtered/", window, "/", window, ".fa.treefile")
        if (!file.exists(fn_tree)) {
            return(data.table::data.table(chr=chr, wsize=window_size, topology="-"))
        }

        # open treefile
        tre <- ape::read.tree(fn_tree)
        tre$edge.length <- NULL
        tre$node.label <- NULL
        
        return(data.table::data.table(chr=chr, wsize=window_size, topology=ape::write.tree(tre)))
    }

    stopCluster(nwcl)

    return(df_output_chr)
}

# function: plot the distribution of window topologies
f_plot_perwindow_topology_dist <- function(df_topology_count, df_topology_count_avg, topology) {
    # create a plot of all window sizes
    df_topology_count_wsize <- df_topology_count[df_topology_count$wsize!="best",]

    plot <- ggplot(df_topology_count_wsize, aes(x=fct_rev(wsize), y=count_percentage, group=chr, ymin=0, ymax=1)) +
        geom_line(aes(alpha=0.2), size=4) +
        geom_point(aes(colour=chr), size=15) +
        viridis::scale_color_viridis(discrete = TRUE) + 
        labs(colour="Chromosomes")

    if (!is.null(df_topology_count_avg)) {
        df_topology_count_avg_wsize <- df_topology_count_avg[df_topology_count_avg$wsize!="best",]

        plot <- plot + geom_point(data=df_topology_count_avg_wsize, aes(x=wsize, y=count_percentage), colour="red", size=15, shape=18) +
            geom_line(data=df_topology_count_avg_wsize, aes(x=wsize, y=count_percentage), colour="red", size=6)
    }
        
    plot <- plot + ggtitle(paste("Distribution of", topology, "across Window Sizes")) + ylab("Proportion") + xlab("Window Size") +
        guides(alpha="none", size="none", colour=guide_legend(override.aes=list(size=5))) +
        theme(
            plot.title = element_text(hjust = 0.5, size = 50),
            plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),
            axis.title.x=element_text(size = 40, margin = margin(t=20, r=0, b=0, l=0)),
            axis.title.y=element_text(size = 40, margin = margin(t=0, r=20, b=0, l=0)),
            axis.text.y=element_text(size=40),
            axis.text.x=element_text(size=40),
            strip.text=element_text(size=40),
            legend.title=element_text(size=30),
            legend.text=element_text(size=30),
            legend.key.size=unit(2,"cm")
        )
    
    return(plot)
}

# function: plot the distribution of best window topologies
f_plot_bestwindow_topology_dist <- function(df_topology_count) {
    plot <- ggplot(df_topology_count, aes(x=topology, y=count_percentage)) +
        geom_boxplot(lwd=3, outlier.shape = NA) +
        geom_point(aes(colour=chr), size=15) +
        viridis::scale_color_viridis(discrete = TRUE) +
        ggtitle("Topology Distribution for Best Window Sizes across Chromosomes") +
        xlab("Topology") + ylab("Proportion") + labs(colour="Chromosomes") +
        guides(size="none", colour=guide_legend(override.aes=list(size=5))) +
        theme(
            plot.title = element_text(hjust = 0.5, size = 50),
            plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),
            axis.title.x=element_text(size = 40, margin = margin(t=20, r=0, b=0, l=0)),
            axis.title.y=element_text(size = 40, margin = margin(t=0, r=20, b=0, l=0)),
            axis.text.y=element_text(size=40),
            axis.text.x=element_text(size=40),
            strip.text = element_text(size=40),
            legend.title=element_text(size=30),
            legend.text=element_text(size=30),
            legend.key.size=unit(2,"cm")
        )

    return(plot)
}