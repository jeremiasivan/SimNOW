library(doSNOW)

#################################
prefix <- "sim"
nthread <- 50

# general
rmddir <- "~/SimNOW/rmd"
outdir <- "~/simulation"
thread <- 10
redo <-  TRUE

# iqtree2 and seqkit
seqkitdir <- "~/seqkit"
iqtree2dir <- "~/iqtree-2.2.2.2-Linux/bin/iqtree2"
masthmmdir <- "~/iqtree-2.2.5.hmmster-Linux/bin/iqtree2"

ic_type <- "AIC"
mast_model <- "JC+T"
#################################

if (nthread / thread < 1) {
  stop("Invalid number of threads (nthread) vs thread per run (thread)")
}

alldirs <- list.dirs(outdir, full.names = F, recursive = F)
dirname <- grepl(paste("^",prefix,"_",sep=""), alldirs)

# create sets of parameters
reports <- list()
for (i in alldirs[dirname]) {
  out <- paste0(outdir,"/",i,"/",i,".wh.html")
  templist <- list(out=out, params=list(prefix=i,
                                        rmddir=rmddir, outdir=outdir, thread=thread, redo=redo,
                                        seqkitdir=seqkitdir, iqtree2dir=iqtree2dir, masthmmdir=masthmmdir,
                                        ic_type=ic_type,
                                        mast_model=mast_model
  ))
  
  reports <- append(reports, list(templist))
}

# function to create reports for independent run
make_report <- function(r) {
  tf <- tempfile()
  dir.create(tf)
  
  rmarkdown::render(input=paste(rmddir,"/wh_main.Rmd", sep=""),
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
#################################