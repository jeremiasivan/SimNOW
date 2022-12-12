library(doParallel)

#################################
prefix <- "wrap"
nthread <- 2

# general
rmddir <- "~/Documents/SimNOW/rmd"
outdir <- "~/Documents/simulation/rmd_test"

# mast and hmm
iqtree2dir <- "~/Downloads/iqtree-2.2.2.2-MacOSX/bin/iqtree2"
mast_model <- c("JC+FO+T", "JC+T")
mast_tops <- 1

#################################

alldirs <- list.dirs(outdir, full.names = F, recursive = F)
dirname <- grepl(paste("^",prefix,"_",sep=""), alldirs)

# create sets of parameters
reports <- list()
for (i in alldirs[dirname]) {
  for (j in mast_model) {
    out <- paste(outdir,"/",i,"/",gsub("[[:punct:]]","",j),".masthmm.html", sep="")
    templist <- list(out=out, params=list(prefix=i,
                                          rmddir=rmddir, outdir=outdir,
                                          iqtree2dir=iqtree2dir,
                                          mast_model=j, mast_tops=mast_tops
                                          ))
  
    reports <- append(reports, list(templist))
  }
}

# function to create reports for independent run
make_report <- function(r) {
  tf <- tempfile()
  dir.create(tf)
  
  rmarkdown::render(input=paste(rmddir,"/masthmm.Rmd", sep=""),
                    output_file=r$out,
                    intermediates_dir=tf,
                    params=r$params,
                    quiet=TRUE)
  unlink(tf)
}

cl <- parallel::makeCluster(nthread)
doParallel::registerDoParallel(cl)

foreach(r=reports, .errorhandling = 'pass') %dopar% make_report(r)

parallel::stopCluster(cl)