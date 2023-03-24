library(doSNOW)

#################################
ms_r <- c(1.7,17,170,1700)
nreps <- 10
prefix <- "hl"
nthread <- 10

# general
rmddir <- "~/Documents/SimNOW/rmd"
outdir <- "~/Documents/simulation/hl"
thread <- 5
redo <- FALSE
  
# sequence simulation
msdir <- "~/Documents/msdir/ms"
ms_params <- "7 1 -T -I 7 1 1 1 1 1 1 1 -ej 2.13 2 1 -ej 5.60 3 1 -ej 7.07 4 1 -ej 6.67 5 6 -ej 8.13 6 1 -ej 15.47 7 1 -es 0.21 2 0.966 -ej 0.21 8 1 -es 0.21 1 0.796 -ej 0.21 9 2 -es 5.52 3 0.115 -ej 5.52 10 4 -es 6.84 6 0.778 -ej 6.84 11 4 -es 6.84 4 0.853 -ej 6.84 12 6"
ms_l <- 10000000
  
iqtree2dir <- "~/Downloads/iqtree-2.2.2.2-MacOSX/bin/iqtree2"
alisim_model <- "JC"
alisim_scale <- 0.005
outgroup <- "7"
  
# non-overlapping window analysis
window_size <- c(1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000)

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