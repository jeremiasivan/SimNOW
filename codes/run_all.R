library(doSNOW)

#################################
ms_r <- c(0,20,200,2000)
nreps <- 10
prefix <- "sim"
nthread <- 50

# general
codedir <- "~/SimNOW/codes"
outdir <- "~/simulation"
thread <- 10
redo <- FALSE

# sequence simulation
msdir <- "~/msdir/ms"
ms_params <- "7 1 -T -I 7 1 1 1 1 1 1 1 -ej 2.13 2 1 -ej 5.60 3 1 -ej 7.07 4 1 -ej 6.67 5 6 -ej 8.13 6 1 -ej 15.47 7 1 -es 0.21 2 0.966 -ej 0.21 8 1 -es 0.21 1 0.796 -ej 0.21 9 2 -es 5.52 3 0.115 -ej 5.52 10 4 -es 6.84 6 0.778 -ej 6.84 11 4 -es 6.84 4 0.853 -ej 6.84 12 6"
ms_l <- 10000000

iqtree2dir <- "~/iqtree2"
alisim_model <- "JC"
alisim_scale <- 0.005

copy_gaps <- TRUE
src_aln <- "~/empirical_aln.fa"

# non-overlapping window analysis
set_model <- FALSE
set_blmin <- TRUE

dna_model <- alisim_model
outgroup <- "7"

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
  prex <- paste0(prefix,"_",temp_table$id[i])
  out <- paste0(outdir,"/",prex,"/",prex,".html")
  
  currentdir <- paste(outdir,"/",prex, sep="")
  if (!file.exists(currentdir)) {
    dir.create(currentdir, recursive = T)
  }
  
  tempsim <- list(out=out, params=list(prefix=prex,
                                       codedir=codedir, outdir=outdir, redo=redo,
                                       msdir=msdir, ms_params=ms_params, ms_r=temp_table$rrate[i], ms_l=ms_l,
                                       iqtree2dir=iqtree2dir, alisim_model=alisim_model, alisim_scale=alisim_scale,
                                       copy_gaps=copy_gaps, src_aln=src_aln
                                       ))
  
  tempnow <- list(out=out, params=list(prefix=prex,
                                       codedir=codedir, outdir=outdir, thread=thread, redo=redo,
                                       ms_l=ms_l,
                                       iqtree2dir=iqtree2dir, set_model=set_model, set_blmin=set_blmin, dna_model=dna_model, outgroup=outgroup,
                                       window_size=window_size
                                       ))

  repsim <- append(repsim, list(tempsim))
  repnow <- append(repnow, list(tempnow))
}

# function to create reports for independent run
make_repsim <- function(r) {
  tf <- tempfile()
  dir.create(tf)
  
  rmarkdown::render(input=paste0(codedir,"/1_sequence_simulation/1_main.Rmd"),
                    output_file=r$out,
                    intermediates_dir=tf,
                    params=r$params,
                    quiet=TRUE)
  unlink(tf)
}

make_repnow <- function(r) {
  tf <- tempfile()
  dir.create(tf)
  
  rmarkdown::render(input=paste0(codedir,"/2_non_overlapping_window/1_main.Rmd"),
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
rmarkdown::render(input=paste(codedir,"/3_all_runs_summary/1_main.Rmd", sep=""),
                  output_file=paste(outdir, "/", prefix, ".html", sep=""),
                  params=list(prefix=prefix, codedir=codedir, outdir=outdir),
                  quiet=TRUE)

#################################