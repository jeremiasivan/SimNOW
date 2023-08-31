#################################
nthread <- 50

# general
codedir <- "~/SimNOW/codes"
thread <- 10
outdir <- ""
redo <- FALSE

# data_edelman.Rmd
fn_hal <- ""
fn_refseq <- ""

dir_hal2maf <- "~/hal2maf"
dir_singleCopy <- "~/getSingleCopy.py"
dir_mafsort <- "~/maf-sort.sh"
dir_msaview <- "~/msa_view"

# nw_main.Rmd
prefix <- ""
dir_iqtree2 <- "~/iqtree2"

set_blmin <- TRUE
set_model <- TRUE
dna_model <- "JC"
outgroup <- "HmelRef"

window_len <- 20
window_size <- c(50000)
min_window_size <- 1000

window_size_nogaps <- c(25000)
min_window_size_nogaps <- 100

ic_type <- "aic"
#################################

if (nthread / thread < 1) {
  stop("Invalid number of threads (nthread) vs thread per run (thread)")
}

# data conversion from HAL -> FASTA
rmarkdown::render(input=paste0(codedir,"/5_empirical_data/edelman_etal_2019/data_preparation.Rmd"),
                  output_file=paste0(outdir,"/edelman.html"),
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
  out <- paste0(run_outdir,"/",c,"/",c,".html")
  
  currentdir <- paste0(run_outdir,"/",c,"/")
  if (!dir.exists(currentdir)) {
    dir.create(currentdir, recursive = T)
  }
  
  # input_aln <- ifelse(c=="chr_all",
  #                     paste0(outdir,"/",c,"/all_concat_filtered.fa"),
  #                     paste0(outdir,"/",c,"/fasta/concatenation/",c,"_concat_filtered.fa"))

  input_aln <- paste0(outdir,"/",c,"/fasta/concatenation/",c,"_concat_filtered.fa")
  temprun <- list(out=out, params=list(codedir=codedir,
                                       prefix=c, outdir=run_outdir, thread=thread, redo=redo,
                                       iqtree2dir=dir_iqtree2, set_blmin=set_blmin, set_model=set_model, dna_model=dna_model, outgroup=outgroup,
                                       input_aln=input_aln, window_len=window_len, window_size=window_size, min_window_size=min_window_size,
                                       ic_type=ic_type
  )

  input_aln <- paste0(outdir,"/",c,"/fasta/concatenation/",c,"_concat_nogaps.fa")
  temprun <- list(out=out, params=list(codedir=codedir,
                                       prefix=paste0("nogaps_",c), outdir=run_outdir, thread=thread, redo=redo,
                                       iqtree2dir=dir_iqtree2, set_blmin=set_blmin, set_model=set_model, dna_model=dna_model, outgroup=outgroup,
                                       input_aln=input_aln, window_len=window_len, window_size=window_size_nogaps, min_window_size=min_window_size_nogaps,
                                       ic_type=ic_type                                     
  ))
  
  runs <- append(runs, list(temprun))
}

# function to create reports for independent run
make_runs <- function(r) {
  tf <- tempfile()
  dir.create(tf)
  
  rmarkdown::render(input=paste0(codedir,"/5_empirical_data/1_main.Rmd"),
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
rmarkdown::render(input=paste0(codedir,"/5_empirical_data/edelman_etal_2019/summary_all.Rmd"),
                  output_file=paste0(outdir,"/",prefix,"/edelman_summary.html"),
                  params=list(prefix=prefix, outdir=outdir, ic_type=ic_type),
                  quiet=TRUE)

#################################