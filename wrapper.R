library(doParallel)

#################################
# wrapper-specific parameters
ms_r <- c(0,1,10)
nreps <- 2
prefix <- "wrap"
nthread <- 2

# general
rmddir <- "~/Documents/SimNOW/rmd"
outdir <- "~/Documents/simulation/rmd_test"
  
# sequence simulation
msdir <- "~/Documents/msdir/ms"
ms_params <- "7 1 -T -I 7 1 1 1 1 1 1 1 -ej 10 2 1 -ej 24 3 1 -ej 35 4 1 -ej 32 5 6 -ej 53 6 1 -ej 69 7 1 -es 1 2 0.25 -ej 1 8 1 -es 23 3 0.25 -ej 23 9 4 -es 34 6 0.25 -ej 34 10 4"
ms_l <- 1000000
  
iqtree2dir <- "~/Downloads/iqtree-2.2.2.2-MacOSX/bin/iqtree2"
alisim_model <- "JC"
alisim_scale <- 0.00029
outgroup <- "7"
  
# non-overlapping window analysis
window_size <- c(50000,100000)

#################################

temp_table <- expand.grid(seq = seq(1:nreps), rrate = ms_r)
temp_table$id <- rownames(temp_table)

cl <- parallel::makeCluster(nthread)
doParallel::registerDoParallel(cl)

foreach(i = 1:nrow(temp_table), .errorhandling = 'pass') %dopar% {
  prex <- paste(prefix,"_",temp_table$id[i], sep="")
  rmdmaindir <- paste(rmddir,"/main.Rmd", sep="")
  
  rmarkdown::render(rmdmaindir, params=list(prefix=prex,
                                            rmddir=rmddir, outdir=outdir,
                                            msdir=msdir, ms_params=ms_params, ms_r=temp_table$rrate[i], ms_l=ms_l,
                                            iqtree2dir=iqtree2dir, alisim_model=alisim_model, alisim_scale=alisim_scale, outgroup=outgroup,
                                            window_size=window_size
                                            ))
}

parallel::stopCluster(cl)