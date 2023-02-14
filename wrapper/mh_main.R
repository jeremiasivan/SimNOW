library(doSNOW)

#################################
prefix <- "hl"
nthread <- 4

# general
rmddir <- "~/Documents/SimNOW/rmd"
outdir <- "~/Documents/simulation/hl"
redo <-  TRUE

# mast and hmm
iqtree2dir <- "~/Downloads/iqtree-2.2.2.2-MacOSX/bin/iqtree2"
analysis_type <- "HMMSTER"
mast_model <- c("JC+T","JC+G+T")
mast_tops <- 0.95

#################################

alldirs <- list.dirs(outdir, full.names = F, recursive = F)
dirname <- grepl(paste("^",prefix,"_",sep=""), alldirs)

# create sets of parameters
reports <- list()
for (i in alldirs[dirname]) {
  for (j in mast_model) {
    out <- paste(outdir,"/",i,"/",gsub("[[:punct:]]","",j),".",tolower(analysis_type),".html", sep="")
    templist <- list(out=out, params=list(prefix=i,
                                          rmddir=rmddir, outdir=outdir, redo=redo,
                                          iqtree2dir=iqtree2dir,
                                          analysis_type=analysis_type, mast_model=j, mast_tops=mast_tops
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
cl <- makeCluster(nthread)
registerDoSNOW(cl)

foreach(r=reports, .errorhandling = 'pass') %dopar% make_report(r)

stopCluster(cl)

# summary
rmarkdown::render(input=paste(rmddir,"/summary_mh.Rmd", sep=""),
                  output_file=paste(outdir, "/", prefix, ".mh.html", sep=""),
                  params=list(prefix=prefix, outdir=outdir, redo=redo, analysis_type=analysis_type),
                  quiet=TRUE)

#################################