# SimNOW

## Table of Content
- <a href="#genpipe">General Pipeline</a>
- <a href="#foldstr">Folder Structure</a>
- <a href="#analyses">Analyses</a>
    - <a href="#seqsim">Sequence simulation</a>
    - <a href="#now">Non-overlapping window analysis</a>
    - <a href="#hmm">HMM</a>
- <a href="#notes">Notes</a>
    - <a href="#example">Example from publication</a>

## <a id="genpipe">General Pipeline</a>
1. **Clone the Git repository** <br>
    ```
    git clone git@github.com:jeremiasivan/SimNOW.git
    ```
2. **Install the prerequisites** <br>
    All software have to be executable, while the R packages should be installed either in your local directory or virtual environment. We recommend you to use environment management system (e.g. `conda`) to install the packages, but you can also use `install.packages()` built-in function in R or RStudio.

    For setting up `conda` environment, please visit the <a href="https://conda.io/projects/conda/en/latest/user-guide/index.html">user guide</a>.
3. **Update the parameters in the code file** <br>
    The details of each parameter is described below.
4. **Run the code file** <br>
    The details for running the code is described below.

## <a id="foldstr">Folder Structure</a>
The following chart shows the general folder structure after running the first two steps. Explanation of each output file is described below.
```
prefix/
├── simulation/
│   ├── prefix.ms
│   ├── prefix.mstrees
│   ├── prefix.mstops
│   ├── prefix.uqtops
│   ├── prefix.nex
│   ├── prefix.multitree
│   ├── prefix.fa
│   └── prefix_nogaps.fa (optional)
├── windows/
│   ├── 10000/
│   │   ├── alignment/
│   │   │   ├── window_01.fa
│   │   │   ├── window_02.fa
│   │   │   ├── window_03.fa
│   │   │   ├── window_04.fa
│   │   │   └── window_05.fa
│   │   ├── summary/
│   │   │   ├── prefix.atsum
│   │   │   ├── prefix.cmp
│   │   │   ├── prefix.cmptw
│   │   │   ├── prefix.cnsum
│   │   │   └── prefix.topsum
│   │   └── trees/
│   │       ├── prefix.best_model.nex
│   │       ├── prefix.ckp.gz
│   │       ├── prefix.iqtree
│   │       ├── prefix.log
│   │       ├── prefix.parstree
│   │       └── prefix.treefile
│   ├── 25000/
│   │   └── ...
│   ├── 50000/
│   │   └── ...
│   ├── prefix.aic.tiff
│   ├── prefix.bic.tiff
│   ├── prefix.aicc.tiff
│   └── prefix.sum
├── prefix.html
└── prefix.log 
```

## <a id="analyses">Analyses</a>
> **Important Notes** <br>
In this pipeline, the results of one step might be used in the other steps. Please be careful about the *consistency* of the naming convention and folder structure.

### <a id="seqsim">Sequence simulation</a>
In this step, we simulate locus trees using `ms` which then inputted into `AliSim` to generate the simulated alignment. The parameters for this step is set in `1_sequence_simulation/simulation.Rmd`.

| Parameters     | Definition                                                                             |
| -------------- | -------------------------------------------------------------------------------------- |
| `prefix`       | Prefix for output files and folder                                                     | 
| `outdir`       | Output directory                                                                       |
| `redo`         | If `FALSE`, skip analysis if output files exist; if `TRUE`, overwrite previous results |
| `msdir`        | Directory for `ms` executable                                                          |
| `ms_params`    | Parameters for `ms` **without** recombination rate (`-r`)                              |
| `ms_r`         | Recombination rate for `ms`                                                            |
| `ms_l`         | Alignment length                                                                       |
| `iqtree2dir`   | Directory for `IQ-Tree2` executable with `AliSim`                                      |
| `alisim_model` | DNA substitution model for simulating alignment in `AliSim`                            |
| `alisim_scale` | Branch length scaling factor for `AliSim`                                              |
| `copy_gaps`    | If `TRUE`, copy gap pattern from reference alignment                                   |
| `src_aln`      | Filepath of the reference alignment (evaluated when `copy_gaps` == `TRUE`)             |

#### Output
Running the code will create a new folder called `simulation` with the following files:
- `prefix.ms`: output `ms` file
- `prefix.mstrees`: table of locus boundaries and trees
- `prefix.mstops`: table of locus boundaries and topologies
- `prefix.uqtops`: list of unique topology
- `prefix.nex`: partition file based on `ms` output trees
- `prefix.multitree`: list of locus trees from `prefix.mstrees`
- `prefix.fa`: FASTA alignment from `AliSim`
- `prefix_nogaps.fa`: FASTA alignment without gaps (only when `only_gaps` == `TRUE`)

> **Important Notes** <br>
The term `topology` refers to phylogenetic tree without branch lengths (i.e., only record the speciation events and not the timing).

### <a id="now">Non-overlapping window analysis</a>
In this step, we use the simulated alignment to perform non-overlapping window analysis and generate the summary statistics. The parameters for this step is set in `2_non_overlapping_window/1_main.Rmd`.

| Parameters     | Definition                                                                                                                            |
| -------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `codedir`      | Directory for folder `SimNOW/codes`                                                                                                   |
| `prefix`       | Prefix for output files and folder                                                                                                    | 
| `outdir`       | Output directory                                                                                                                      |
| `thread`       | Number of threads for parallelization                                                                                                 |
| `redo`         | If `FALSE`, skip analysis if output files exist; if `TRUE`, overwrite previous results                                                |
| `ms_l`         | Alignment length                                                                                                                      |
| `iqtree2dir`   | Directory for `IQ-Tree2` executable                                                                                                   |
| `set_model`    | If `TRUE`, set the substitution model of window trees to be `dna_model`; if `FALSE` use ModelFinder to find the best model per window |
| `set_blmin`    | If `TRUE`, set the minimum branch length to be 1/window_size                                                                          |
| `dna_model`    | DNA substitution model for window trees                                                                                               |
| `outgroup`     | Outgroup of window trees (optional)                                                                                                   |
| `window_size`  | Vector of window sizes; each has to be the factor of `ms_l`                                                                           |

#### Output
Running the code will create a new folder called `windows` with individual folder for each window size. Each window size folder consists of three folders:
- `alignment/`: folder with all window alignments
- `summary/`
    - prefix.atsum: table of window boundaries and trees
    - prefix.cmp: matrix of sites that recover specific topology according to `ms` (row) and non-overlapping windows (column)
    - prefix.cmptw: comparison table of `ms` and non-overlapping windows topology distribution
    - prefix.cnsum: frequency of consecutive windows that recover the same topology
    - prefix.topsum: table of unique topologies sorted from the most common one
- `trees/`: folder with window trees generated using `IQ-Tree2`

Additionally, the `windows` folder will contain the following files:
- `prefix.aic.tiff`: correlation plot between AIC and window size
- `prefix.bic.tiff`: correlation plot between BIC and window size
- `prefix.aicc.tiff`: correlation plot between AICc and window size
- `prefix.sum`: summary table of window sizes and their respective information criteria scores

### <a id="hmm">HMM</a>
In this step, we use MAST-HMM to refine the breakpoints between windows. We simplify the process by running MAST-HMM on two consecutive windows if the windows have different topology. The parameters for this step is set in `3_hmm/1_main.Rmd`.

| Parameters   | Definition                                                                             |
| ------------ | -------------------------------------------------------------------------------------- |
| `codedir`    | Directory for folder `SimNOW/codes`                                                    |
| `prefix`     | Prefix for output files and folder                                                     | 
| `outdir`     | Output directory                                                                       |
| `thread`     | Number of threads for parallelization                                                  |
| `redo`       | If `FALSE`, skip analysis if output files exist; if `TRUE`, overwrite previous results |
| `seqkitdir`  | Directory for `Seqkit` executable                                                      |
| `iqtree2dir` | Directory for `IQ-Tree2` executable                                                    |
| `masthmmdir` | Directory for `IQ-Tree2` executable with `MAST-HMM`                                    |
| `ic_type`    | Information criterion (i.e., `AIC`, `BIC`, or `AICc`)                                  |
| `mast_model` | DNA substitution model for `MAST` analysis                                             |

#### Output
Running the code will create a new folder called `winhmm` with `MAST-HMM` output. *Further details will be provided in the near future.*

## <a id="notes">Notes</a>
If you want to run more than one analysis at the same time, you can run `run_all.R`. There are several parameters that are new or changed:
| Parameters | Definition                                                                       |
| ---------- | -------------------------------------------------------------------------------- |
| `ms_r`     | Vector of `ms` recombination rate                                                |
| `nreps`    | Number of replicate per `ms` recombination rate                                  |
| `nthread`  | Total number of threads; it is divided by `thread` to parallelize the replicates |

### <a id="example">Example from publication</a>
In the paper, we ran three different simulation scenarios with different degree of incomplete lineage sorting (ILS) but the same percentage of informative sites. Please do refer to the publication for more detailed explanation.

Parameters (in `codes/run_all.R`) with consistent value across scenarios:
```
nreps <- 10
nthread <- 50

codedir <- "~/SimNOW/codes"
thread <- 10
redo <- FALSE

msdir <- "~/msdir/ms"
ms_l <- 10000000

iqtree2dir <- "~/iqtree-2.2.2.4-Linux/bin/iqtree2"
alisim_model <- "JC"
copy_gaps <- FALSE
src_aln <- ""

set_model <- TRUE
set_blmin <- FALSE
dna_model <- alisim_model
outgroup <- "7"

window_size <- c(100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000)

seqkitdir <- "~/seqkit"
masthmmdir <- "~/iqtree-2.2.5.hmmster-Linux/bin/iqtree2"
mast_model <- "JC+T"

ic_type <- "AIC"
```

Specific parameters (in `codes/run_all.R`) for each scenario: <br>
*Low ILS*
```
ms_r <- c(0,2,20,200)
prefix <- "low"
outdir <- "~/Documents/low"

ms_params <- "7 1 -T -I 7 1 1 1 1 1 1 1 -ej 16.0 2 1 -ej 42.0 3 1 -ej 53.0 4 1 -ej 50.0 5 6 -ej 61.0 6 1 -ej 116.0 7 1 -es 1.6 2 0.966 -ej 1.6 8 1 -es 1.6 1 0.796 -ej 1.6 9 2 -es 41.4 3 0.115 -ej 41.4 10 4 -es 51.3 6 0.778 -ej 51.3 11 4 -es 51.3 4 0.853 -ej 51.3 12 6"
alisim_scale <- 0.0007
```

*Mid ILS*
```
ms_r <- c(0,20,200,2000)
prefix <- "mid"
outdir <- "~/Documents/mid"

ms_params <- "7 1 -T -I 7 1 1 1 1 1 1 1 -ej 2.13 2 1 -ej 5.60 3 1 -ej 7.07 4 1 -ej 6.67 5 6 -ej 8.13 6 1 -ej 15.47 7 1 -es 0.21 2 0.966 -ej 0.21 8 1 -es 0.21 1 0.796 -ej 0.21 9 2 -es 5.52 3 0.115 -ej 5.52 10 4 -es 6.84 6 0.778 -ej 6.84 11 4 -es 6.84 4 0.853 -ej 6.84 12 6"
alisim_scale <- 0.005
```

*High ILS*
```
ms_r <- c(0,30,300,3000)
prefix <- "high"
outdir <- "~/Documents/high"

ms_params <- "7 1 -T -I 7 1 1 1 1 1 1 1 -ej 0.16 2 1 -ej 0.42 3 1 -ej 0.53 4 1 -ej 0.50 5 6 -ej 0.61 6 1 -ej 1.16 7 1 -es 0.016 2 0.966 -ej 0.016 8 1 -es 0.016 1 0.796 -ej 0.016 9 2 -es 0.414 3 0.115 -ej 0.414 10 4 -es 0.513 6 0.778 -ej 0.513 11 4 -es 0.513 4 0.853 -ej 0.513 12 6"
alisim_scale <- 0.041
```

---
*Last update: 30 July 2023 by Jeremias Ivan*