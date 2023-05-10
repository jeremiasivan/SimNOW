library(doSNOW)

#################################
prefix <- "sim"
nthread <- 50

# general
rmddir <- "~/SimNOW/rmd"
outdir <- "~/simulation"
thread <- 10
redo <-  TRUE

# mast and hmm
iqtree2dir <- "~/iqtree-2.2.2.4-MacOSX/bin/iqtree2"

ic_type <- "AIC"

mast_analysis <- "HiMAST"
mast_model <- c("JC+T")
percent_tops <- 0.95
max_tops <- 10
#################################

if (nthread / thread < 1) {
  stop("Invalid number of threads (nthread) vs thread per run (thread)")
}

alldirs <- list.dirs(outdir, full.names = F, recursive = F)
dirname <- grepl(paste("^",prefix,"_",sep=""), alldirs)

# create sets of parameters
reports <- list()
for (i in alldirs[dirname]) {
  for (j in mast_model) {
    out <- paste(outdir,"/",i,"/",tolower(gsub("[[:punct:]]","",j)),".",tolower(mast_analysis),".html", sep="")
    templist <- list(out=out, params=list(prefix=i,
                                          rmddir=rmddir, outdir=outdir, thread=thread, redo=redo,
                                          iqtree2dir=iqtree2dir,
                                          ic_type=ic_type,
                                          mast_analysis=mast_analysis, mast_model=j, percent_tops=percent_tops, max_tops=max_tops
                                          ))
  
    reports <- append(reports, list(templist))
  }
}

# function to create reports for independent run
make_report <- function(r) {
  tf <- tempfile()
  dir.create(tf)
  
  rmarkdown::render(input=paste(rmddir,"/mh_main.Rmd", sep=""),
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
rmarkdown::render(input=paste(rmddir,"/summary_mh.Rmd", sep=""),
                  output_file=paste(outdir, "/", prefix, ".mh.html", sep=""),
                  params=list(prefix=prefix, outdir=outdir, redo=redo, mast_analysis=mast_analysis),
                  quiet=TRUE)

#################################