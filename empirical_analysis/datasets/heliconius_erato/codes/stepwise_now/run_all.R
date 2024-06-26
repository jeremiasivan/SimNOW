library(doSNOW)

#################################
nthread <- 50

# general
codedir <- "~/SimNOW/empirical_analysis"
thread <- 10
outdir <- ""
redo <- FALSE

# data_preparation.Rmd
fn_hal <- ""
fn_refseq <- paste0(codedir, "/datasets/heliconius_erato/files/refseq.txt")

dir_hal2maf <- "~/hal2maf"
dir_singleCopy <- "~/getSingleCopy.py"
dir_mafsort <- "~/maf-sort.sh"
dir_msaview <- "~/msa_view"

# stepwise_now
prefix <- ""
dir_iqtree2 <- "~/iqtree2"

set_blmin <- FALSE
set_model <- FALSE
dna_model <- ""

bootstrap <- 1000
bootstrap_type <- "ufboot"
outgroup <- "HmelRef"

initial_wsize <- 64000
min_wsize <- 0
min_informative_sites <- 4

# summary
colour_scheme <- paste0(codedir, "/datasets/heliconius_erato/files/colour_scheme.txt")

#################################

if (nthread / thread < 1) {
  stop("Invalid number of threads (nthread) vs thread per run (thread)")
}

# data conversion from HAL -> FASTA
rmarkdown::render(input=paste0(codedir,"/datasets/heliconius_erato/codes/data_preparation.Rmd"),
                  output_file=paste0(outdir,"/heliconius.html"),
                  params=list(fn_hal=fn_hal, fn_refseq=fn_refseq, thread=thread, outdir=outdir, redo=redo,
                              dir_hal2maf=dir_hal2maf, dir_singleCopy=dir_singleCopy, dir_mafsort=dir_mafsort, dir_msaview=dir_msaview),
                  quiet=TRUE)

# run NOW on every chromosome
ls_chr <- list.dirs(outdir, full.names=F, recursive=F)
ls_chr <- subset(ls_chr, grepl("chr+", ls_chr))

# create sets of parameters
run_outdir <- paste0(outdir,"/",prefix,"/")
runs <- list()

for (c in ls_chr) {
  currentdir <- paste0(run_outdir,"/",c,"/")
  if (!dir.exists(currentdir)) {
    dir.create(currentdir, recursive = T)
  }

  # filtered sequence
  out <- paste0(currentdir,c,".html")
  input_aln <- paste0(outdir,"/",c,"/fasta/concatenation/",c,"_concat.fa")
  temprun <- list(out=out, params=list(codedir=codedir,
                                       prefix=c, outdir=run_outdir, thread=thread, redo=redo,
                                       iqtree2dir=dir_iqtree2, set_blmin=set_blmin, set_model=set_model, dna_model=dna_model, 
                                       bootstrap=bootstrap, bootstrap_type=bootstrap_type, outgroup=outgroup,
                                       input_aln=input_aln, initial_wsize=initial_wsize, min_wsize=min_wsize, min_informative_sites=min_informative_sites
  ))
  runs <- append(runs, list(temprun))
}

# function to create reports for independent run
make_runs <- function(r) {
  tf <- tempfile()
  dir.create(tf)
  
  rmarkdown::render(input=paste0(codedir,"/stepwise_now/1_main.Rmd"),
                    output_file=r$out,
                    intermediates_dir=tf,
                    params=r$params,
                    quiet=TRUE)
  unlink(tf)
}

# run parallelized EmpNOW analysis
cl <- makeCluster(floor(nthread/thread), outfile="")
registerDoSNOW(cl)

foreach(r=runs, .errorhandling = 'pass') %dopar% {
  make_runs(r)
  NULL
}

stopCluster(cl)

# summary for all chromosomes
rmarkdown::render(input=paste0(codedir,"/datasets/heliconius_erato/codes/stepwise_now/summary_all.Rmd"),
                  output_file=paste0(outdir,"/",prefix,"/heliconius_summary.html"),
                  params=list(codedir=codedir, prefix=prefix, thread=thread, outdir=outdir, initial_wsize=initial_wsize, colour_scheme=colour_scheme),
                  quiet=TRUE)
                  
#################################