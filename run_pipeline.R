#!/usr/bin/env Rscript

# ============================================================
#  SimNOW
#
#  Usage: Rscript run_pipeline.R --config config.yaml
#         Rscript run_pipeline.R --config config.yaml --redo
# ============================================================

# --- Load libraries and function ----------------------------
suppressPackageStartupMessages(library(doSNOW))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(yaml))

# create a function to retrieve parameter value with default
f_get_param <- function(value, default) {
  if (is.null(value) || identical(value, "")) {
    default
  } else {
    value
  }
}

# --- Argument parsing ----------------------------------------
option_list <- list(
  make_option(c("-c", "--config"), type="character", default=NULL,
              help="Path to YAML config file [required]", metavar="FILE"),
  make_option(c("-r", "--redo"), action="store_true", default=FALSE,
              help="Re-run all analyses and override previous results")
)

# parse the arguments
opt <- parse_args(OptionParser(option_list=option_list))

# stop the code if config file is invalid
if (is.null(opt$config)) {
  stop("-c/--config is required.\n
        Usage: Rscript run_pipeline.R --config config.yaml")
}

if (!file.exists(opt$config)) {
  stop(paste("Config file not found:", opt$config))
}

# --- Load and validate config --------------------------------
cfg <- yaml::read_yaml(opt$config)

# set required parameters
required_fields <- c("codedir", "outdir", "exe_iqtree", "exe_ms", "exe_seqkit", "fn_ms_metadata", "ms_params", "window_size")
missing <- setdiff(required_fields, names(cfg))
if (length(missing) > 0) {
  stop(paste("Missing required config fields:", paste(missing, collapse=", ")))
}

# check if input ms recombination rates file is invalid
if (is.null(cfg$fn_ms_metadata) || cfg$fn_ms_metadata == "" || !file.exists(cfg$fn_ms_metadata)) {
  stop(paste("fn_ms_r file not found:", cfg$fn_ms_r))
}

# check if window_size is more than one
if (length(cfg$window_size) <= 1) {
  stop(paste("Invalid set of window sizes:", cfg$window_size))
}

# apply defaults for substitution models
cfg$alisim_model <- f_get_param(cfg$alisim_model, "JC")
cfg$dna_model <- f_get_param(cfg$dna_model, cfg$alisim_model)
cfg$nreps <- as.integer(f_get_param(cfg$nreps, 10))

# --- Apply CLI overrides -------------------------------------
if (opt$redo) {
  message("Note: --redo flag set via CLI, overriding config.")
  cfg$redo <- TRUE
}

# --- Map config to rmarkdown params --------------------------
render_params <- list(
  codedir               = cfg$codedir,
  prefix                = f_get_param(cfg$prefix, "simulation"),
  outdir                = cfg$outdir,
  thread                = as.integer(f_get_param(cfg$thread, 1)),
  redo                  = as.logical(f_get_param(cfg$redo, FALSE)),

  exe_iqtree            = cfg$exe_iqtree,
  exe_ms                = cfg$exe_ms,
  exe_seqkit            = cfg$exe_seqkit,

  ms_params             = cfg$ms_params,

  alisim_model          = cfg$alisim_model,
  alisim_scale          = as.numeric(f_get_param(cfg$alisim_scale, 0.01)),
  is_copy_gaps          = as.logical(f_get_param(cfg$is_copy_gaps, FALSE)),
  src_aln               = ifelse(file.exists(cfg$src_aln), cfg$src_aln, ""),

  set_blmin             = as.logical(f_get_param(cfg$set_blmin, TRUE)),
  set_model             = as.logical(f_get_param(cfg$set_model, FALSE)),
  set_keepident         = as.logical(f_get_param(cfg$set_keepident, FALSE)),
  dna_model             = cfg$dna_model,
  outgroup              = f_get_param(cfg$outgroup, ""),

  window_size           = cfg$window_size,
  split_prop            = unlist(f_get_param(cfg$split_prop, list(0.25, 0.5, 0.75))),
  min_informative_sites = as.integer(f_get_param(cfg$min_informative_sites, 1))
)

# --- Map config to rmarkdown params --------------------------

# open ms metadata file
df_ms <- data.table::fread(cfg$fn_ms_metadata)
ls_group <- unique(df_ms$group)

# iterate over groupings
temp_table <- data.frame()
for (group in ls_group) {
  df_ms_subset <- df_ms[df_ms$group==group, ]

  temp_table <- rbind(temp_table,
                      data.frame(group  = group,
                                 seq    = 1:cfg$nreps,
                                 rrate  = I(replicate(cfg$nreps, df_ms_subset$ms_r, simplify=FALSE)),
                                 seqlen = I(replicate(cfg$nreps, df_ms_subset$ms_l, simplify=FALSE))))
}

# add ID
temp_table$id <- rownames(temp_table)

# build per-replicate parameter list
repsim <- vector("list", nrow(temp_table))
repfix <- vector("list", nrow(temp_table))
repvar <- vector("list", nrow(temp_table))

for (i in 1:nrow(temp_table)) {
  prex <- paste0(render_params$prefix, "_", temp_table$id[i])
  currentdir <- file.path(render_params$outdir, prex)
  out <- file.path(currentdir, paste0(prex, ".html"))
  
  if (!dir.exists(currentdir)) {
    dir.create(currentdir, recursive = T)
  }

  # extract alignment length
  init_wsize_concat <- sum(unlist(temp_table$seqlen[[i]]))
  init_wsize_concat <- f_get_param(cfg$init_wsize, init_wsize_concat)
  
  repsim[[i]] <- list(out=out,
                      params=c(render_params[c("codedir", "outdir", "redo", "exe_ms", "ms_params",
                                               "exe_iqtree", "alisim_model", "alisim_scale", "is_copy_gaps", "src_aln")],
                               list(prefix=prex, ms_r=unlist(temp_table$rrate[[i]]), ms_l=unlist(temp_table$seqlen[[i]]))
                      ))

  
  repfix[[i]] <- list(out=out,
                      params=c(render_params[c("codedir", "outdir", "thread", "redo", "exe_iqtree",
                                               "set_model", "set_blmin", "dna_model", "outgroup", "window_size")],
                               list(prefix=prex, ms_l=unlist(temp_table$seqlen[[i]]))
                      ))

  repvar[[i]] <- list(out=out,
                      params=c(render_params[c("codedir", "outdir", "thread", "redo", "exe_seqkit", "exe_iqtree",
                                               "set_model", "set_keepident", "set_blmin", "dna_model", "outgroup",
                                               "division_prop", "min_informative_sites")],
                      list(prefix=prex, init_wsize=init_wsize_concat)
                      ))
}

# function to create reports for independent run
make_report <- function(r, rmd_path) {
  tf <- tempfile()
  dir.create(tf)
  on.exit(unlink(tf, recursive=T))
  
  rmarkdown::render(input=rmd_path,
                    output_file=r$out,
                    intermediates_dir=tf,
                    params=r$params,
                    quiet=TRUE)
}

make_repsim <- function(r) make_report(r, file.path(render_params$codedir, "1_sequence_simulation", "1_main.Rmd"))
make_repfix <- function(r) make_report(r, file.path(render_params$codedir, "2_non_overlapping_window", "1_main.Rmd"))
make_repvar <- function(r) make_report(r, file.path(render_params$codedir, "3_variable_window_size", "1_main.Rmd"))

# run parallelized sequence simulations
cl <- makeCluster(thread, outfile="")
registerDoSNOW(cl)

foreach(r=repsim, .errorhandling='pass') %dopar% {
  make_repsim(r)
  NULL
}

stopCluster(cl)

# run NOW analyses
for (r in repfix) make_repfix(r)
for (r in repvar) make_repvar(r)

# summary
rmarkdown::render(input       = file.path(render_params$codedir, "4_all_runs_summary", "1_main.Rmd"),
                  output_file = file.path(render_params$outdir, paste0(render_params$prefix, ".html")),
                  params      = list(prefix=render_params$prefix, codedir=render_params$codedir, outdir=render_params$outdir),
                  quiet       = TRUE)

#################################