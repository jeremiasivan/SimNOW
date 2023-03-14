library(doSNOW)

#################################
ms_r <- c(0,20,100,200,400,600,800,1000)
nreps <- 5
prefix <- "hl"
nthread <- 5

# general
rmddir <- "~/Documents/SimNOW/rmd"
outdir <- "~/Documents/simulation/debug/hl"
thread <- 1
redo <- FALSE
  
# sequence simulation
msdir <- "~/Documents/msdir/ms"
ms_params <- "7 1 -T -I 7 1 1 1 1 1 1 1 -ej 8.0 2 1 -ej 21.0 3 1 -ej 26.5 4 1 -ej 25.0 5 6 -ej 30.5 6 1 -ej 58.0 7 1 -es 4.0 2 0.5 -ej 4.0 8 1 -es 10.5 3 0.5 -ej 10.5 9 4 -es 25.7 6 0.5 -ej 25.7 10 4"
ms_l <- 10000000
  
iqtree2dir <- "~/Downloads/iqtree-2.2.2.2-MacOSX/bin/iqtree2"
alisim_model <- "JC"
alisim_scale <- 0.0015
outgroup <- "7"
  
# non-overlapping window analysis
window_size <- c(50000)

#################################

if (nthread / thread < 1) {
  stop("Invalid number of threads (nthread) vs thread per run (thread)")
}

temp_table <- expand.grid(seq = seq(1:nreps), rrate = ms_r)
temp_table$id <- rownames(temp_table)

# create sets of parameters
reports <- list()
for (i in 1:nrow(temp_table)) {
  prex <- paste(prefix,"_",temp_table$id[i], sep="")
  out <- paste(outdir,"/",prex,"/",prex,".html", sep="")
  
  currentdir <- paste(outdir,"/",prex, sep="")
  if (!file.exists(currentdir)) {
    dir.create(currentdir, recursive = T)
  }
  
  templist <- list(out=out, params=list(prefix=prex,
                                        rmddir=rmddir, outdir=outdir, thread=thread, redo=redo,
                                        msdir=msdir, ms_params=ms_params, ms_r=temp_table$rrate[i], ms_l=ms_l,
                                        iqtree2dir=iqtree2dir, alisim_model=alisim_model, alisim_scale=alisim_scale, outgroup=outgroup,
                                        window_size=window_size
                                        ))
  
  reports <- append(reports, list(templist))
}

# function to create reports for independent run
make_report <- function(r) {
  tf <- tempfile()
  dir.create(tf)
  
  rmarkdown::render(input=paste(rmddir,"/nw_main.Rmd", sep=""),
                    output_file=r$out,
                    intermediates_dir=tf,
                    params=r$params,
                    quiet=TRUE)
  unlink(tf)
}

# run parallelized simulations
cl <- makeCluster(floor(nthread/thread), outfile="")
registerDoSNOW(cl)

foreach(r=reports, .errorhandling = 'pass') %dopar% {
  make_report(r)
  NULL
}

stopCluster(cl)

# summary
rmarkdown::render(input=paste(rmddir,"/summary_nw.Rmd", sep=""),
                  output_file=paste(outdir, "/", prefix, ".html", sep=""),
                  params=list(prefix=prefix, outdir=outdir),
                  quiet=TRUE)

#################################