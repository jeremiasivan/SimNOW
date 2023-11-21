#################################
# general
codedir <- "~/Documents/SimNOW/empirical_analysis"
prefix <- "primates_chrY"
outdir <- "~/Downloads/mammals"
thread <- 2
redo <- FALSE

# IQ-Tree2
iqtree2dir <- "~/Downloads/iqtree-2.2.2.2-MacOSX/bin/iqtree2"
set_blmin <- FALSE
set_model <- FALSE
dna_model <- ""

bootstrap <- 1000
bootstrap_type <- "ufboot"
outgroup <- ""

# NOW analysis
input_aln <- "~/Downloads/mammals/primates_chrY.fa"
initial_wsize <- 5000
min_wsize <- 1000

min_informative_sites <- 0
#################################

# create outdir
outdir_prefix <- paste0(outdir,"/",prefix,"/")
if (!dir.exists(outdir_prefix)) {
    dir.create(outdir_prefix, recursive=T)
}

# run stepwise NOW
rmarkdown::render(input=paste0(codedir,"/stepwise_now/1_main.Rmd"),
                  output_file=paste0(outdir_prefix,prefix,".html"),
                  params=list(codedir=codedir, prefix=prefix, thread=thread, outdir=outdir, redo=redo,
                              iqtree2dir=iqtree2dir, set_blmin=set_blmin, set_model=set_model, dna_model=dna_model,
                              bootstrap=bootstrap, bootstrap_type=bootstrap_type, outgroup=outgroup,
                              input_aln=input_aln, initial_wsize=initial_wsize, min_wsize=min_wsize, min_informative_sites=min_informative_sites),
                  quiet=TRUE)
