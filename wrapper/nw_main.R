library(doSNOW)

#################################
ms_r <- c(0,20,200,2000)
nreps <- 10
prefix <- "sim"
nthread <- 50

# general
rmddir <- "~/SimNOW/rmd"
outdir <- "~/simulation"
thread <- 10
redo <- FALSE
  
# sequence simulation
msdir <- "~/msdir/ms"
ms_params <- "7 1 -T -I 7 1 1 1 1 1 1 1 -ej 2.13 2 1 -ej 5.60 3 1 -ej 7.07 4 1 -ej 6.67 5 6 -ej 8.13 6 1 -ej 15.47 7 1 -es 0.21 2 0.966 -ej 0.21 8 1 -es 0.21 1 0.796 -ej 0.21 9 2 -es 5.52 3 0.115 -ej 5.52 10 4 -es 6.84 6 0.778 -ej 6.84 11 4 -es 6.84 4 0.853 -ej 6.84 12 6"
ms_l <- 10000000
  
iqtree2dir <- "~/iqtree-2.2.2.4-MacOSX/bin/iqtree2"
alisim_model <- "JC"
alisim_scale <- 0.005
outgroup <- "7"
  
# non-overlapping window analysis
window_size <- c(100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000)

#################################

if (nthread / thread < 1) {
  stop("Invalid number of threads (nthread) vs thread per run (thread)")
}

temp_table <- expand.grid(seq = seq(1:nreps), rrate = ms_r)
temp_table$id <- rownames(temp_table)

# create sets of parameters
repsim <- list()
repnow <- list()

for (i in 1:nrow(temp_table)) {
  prex <- paste(prefix,"_",temp_table$id[i], sep="")
  out <- paste(outdir,"/",prex,"/",prex,".html", sep="")
  
  currentdir <- paste(outdir,"/",prex, sep="")
  if (!file.exists(currentdir)) {
    dir.create(currentdir, recursive = T)
  }
  
  tempsim <- list(out=out, params=list(prefix=prex,
                                       outdir=outdir, redo=redo,
                                       msdir=msdir, ms_params=ms_params, ms_r=temp_table$rrate[i], ms_l=ms_l,
                                       iqtree2dir=iqtree2dir, alisim_model=alisim_model, alisim_scale=alisim_scale
                                       ))
  
  tempnow <- list(out=out, params=list(prefix=prex,
                                       rmddir=rmddir, outdir=outdir, thread=thread, redo=redo,
                                       ms_l=ms_l,
                                       iqtree2dir=iqtree2dir, alisim_model=alisim_model, outgroup=outgroup,
                                       window_size=window_size
                                       ))
  
  repsim <- append(repsim, list(tempsim))
  repnow <- append(repnow, list(tempnow))
}

# function to create reports for independent run
make_repsim <- function(r) {
  tf <- tempfile()
  dir.create(tf)
  
  rmarkdown::render(input=paste(rmddir,"/nw_simulation.Rmd", sep=""),
                    output_file=r$out,
                    intermediates_dir=tf,
                    params=r$params,
                    quiet=TRUE)
  unlink(tf)
}

make_repnow <- function(r) {
  tf <- tempfile()
  dir.create(tf)
  
  rmarkdown::render(input=paste(rmddir,"/nw_main.Rmd", sep=""),
                    output_file=r$out,
                    intermediates_dir=tf,
                    params=r$params,
                    quiet=TRUE)
  unlink(tf)
}

# run parallelized sequence simulations
cl <- makeCluster(nthread, outfile="")
registerDoSNOW(cl)

foreach(r=repsim, .errorhandling = 'pass') %dopar% {
  make_repsim(r)
  NULL
}

stopCluster(cl)

# run parallelized NOW analysis
cl <- makeCluster(floor(nthread/thread), outfile="")
registerDoSNOW(cl)

foreach(r=repnow, .errorhandling = 'pass') %dopar% {
  make_repnow(r)
  NULL
}

stopCluster(cl)

# summary
rmarkdown::render(input=paste(rmddir,"/summary_nw.Rmd", sep=""),
                  output_file=paste(outdir, "/", prefix, ".html", sep=""),
                  params=list(prefix=prefix, outdir=outdir),
                  quiet=TRUE)

#################################