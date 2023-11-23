#################################
nthread <- 50

# general
codedir <- "~/SimNOW/empirical_analysis"
thread <- 10
outdir <- ""
redo <- FALSE

# data_preparation.Rmd
fn_mafseq <- "~/SimNOW/empirical_analysis/datasets/ucsc_primates/files/mafseq.txt"
dir_msaview <- "~/msa_view"

# stepwise_now
prefix <- "primates"
dir_iqtree2 <- "~/iqtree2"

set_blmin <- TRUE
set_model <- FALSE
dna_model <- ""

bootstrap <- 1000
bootstrap_type <- "ufboot"
outgroup <- ""

initial_wsize <- 64000
min_wsize <- 0

min_informative_sites <- 0
#################################

if (nthread / thread < 1) {
  stop("Invalid number of threads (nthread) vs thread per run (thread)")
}

# data conversion from MAF -> FASTA
rmarkdown::render(input=paste0(codedir,"/datasets/ucsc_primates/codes/data_preparation.Rmd"),
                  output_file=paste0(outdir,"/edelman.html"),
                  params=list(fn_mafseq=fn_mafseq, thread=thread, outdir=outdir, redo=redo, dir_msaview=dir_msaview),
                  quiet=TRUE)

# run NOW on every chromosome
ls_chr <- list.files(paste0(outdir,"/data/fasta/primates/"), pattern='*.fa$', full.names=F, recursive=F)
ls_chr <- sapply(ls_chr, function(x) { gsub(".fa", "", x) })

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
  input_aln <- paste0(outdir,"/data/fasta/primates/",c,".fa")
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
                  
#################################