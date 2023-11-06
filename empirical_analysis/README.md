# SimNOW

## Table of Content
- <a href="#analyses">Analyses</a>
    - <a href="#now">Non-overlapping window analysis</a>
    - <a href="#stepnow">Stepwise non-overlapping windows</a>
    - <a href="#hmm">HMM</a>
- <a href="#datasets">Datasets</a>
    - <a href="#edelman">Edelman et al. (2019)</a>

## <a id="analyses">Analyses</a>

### <a id="now">Non-overlapping window analysis</a>
In this step, we run non-overlapping window analysis on empirical alignment and generate the summary statistics. The parameters for this step is set in `1_main.Rmd`.

| Parameters               | Definition                                                                                                                            |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------- |
| `codedir`                | Directory for folder `SimNOW/codes_empirical`                                                                                         |
| `prefix`                 | Prefix for output files and folder                                                                                                    | 
| `outdir`                 | Output directory                                                                                                                      |
| `thread`                 | Number of threads for parallelization                                                                                                 |
| `redo`                   | If `FALSE`, skip analysis if output files exist; if `TRUE`, overwrite previous results                                                |
| `iqtree2dir`             | Directory for `IQ-Tree2` executable                                                                                                   |
| `set_model`              | If `TRUE`, set the substitution model of window trees to be `dna_model`; if `FALSE` use ModelFinder to find the best model per window |
| `set_blmin`              | If `TRUE`, set the minimum branch length to be 1/window_size                                                                          |
| `dna_model`              | DNA substitution model for window trees                                                                                               |
| `outgroup`               | Outgroup of window trees (optional)                                                                                                   |
| `input_aln`              | Input FASTA alignment                                                                                                                 |
| `window_len`             | Number of window sizes                                                                                                                |
| `window_size`            | Vector of window sizes; if any, trailing positions of `input_aln` will be removed to match the least common multiple of the vector    |
| `min_window_size`        | Minimum window size                                                                                                                   |
| `window_size_nogaps`     | Similar with `window_size`, but only applied for `no_gaps` alignments                                                                 |
| `min_window_size_nogaps` | Similar with `min_window_size`, but only applied for `no_gaps` alignments                                                             |

#### Output
Running the code will create a new folder called `windows` with individual folder for each window size. Each window size folder consists of three folders:
- `alignment/`: folder with all window alignments
- `summary/`
    - `prefix.atsum`: table of window boundaries and trees
    - `prefix.cnsum`: frequency of consecutive windows that recover the same topology
    - `prefix.topsum`: table of unique topologies sorted from the most common one
- `trees/`: folder with window trees generated using `IQ-TREE 2`

Additionally, the `windows` folder will contain the following files:
- `prefix.sum`: summary table of window sizes and their respective information criteria scores
- `prefix.aic.tiff`: correlation plot between AIC and window size
- `prefix.bic.tiff`: correlation plot between BIC and window size
- `prefix.aicc.tiff`: correlation plot between AICc and window size
- `prefix.topdist`: summary table of topology distribution per window size
- `prefix.topdist.tiff`: plot of topology distribution per window size

### <a id="stepnow">Stepwise non-overlapping window analysis</a>
In this step, we run stepwise non-overlapping window analysis on empirical alignment and generate the summary statistics. The parameters for this step is set in `1_main.Rmd`.

| Parameters               | Definition                                                                                                                            |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------- |
| `codedir`                | Directory for folder `SimNOW/codes_empirical`                                                                                         |
| `prefix`                 | Prefix for output files and folder                                                                                                    | 
| `outdir`                 | Output directory                                                                                                                      |
| `thread`                 | Number of threads for parallelization                                                                                                 |
| `redo`                   | If `FALSE`, skip analysis if output files exist; if `TRUE`, overwrite previous results                                                |
| `iqtree2dir`             | Directory for `IQ-Tree2` executable                                                                                                   |
| `set_model`              | If `TRUE`, set the substitution model of window trees to be `dna_model`; if `FALSE` use ModelFinder to find the best model per window |
| `set_blmin`              | If `TRUE`, set the minimum branch length to be 1/window_size                                                                          |
| `dna_model`              | DNA substitution model for window trees                                                                                               |
| `bootstrap`              | Number of bootstrap (required if `bootstrap_type != ""`)                                                                              |
| `bootstrap_type`         | Type of boostrap (either `parametric` or `ufboot`) (optional)                                                                         |
| `outgroup`               | Outgroup of window trees (optional)                                                                                                   |
| `input_aln`              | Input FASTA alignment                                                                                                                 |
| `initial_wsize`          | Starting (biggest) window size                                                                                                        |
| `min_informative_sites`  | Minimum number of informative sites per window                                                                                        |

#### Output
Running the code will create a new folder for each pair of window sizes, starting from the `initial_wsize` and its half (i.e., `initial_wsize/2`). Each folder is named incrementally from `1` (hitherto referred as `step`) and comprised of the following folders and files:
- `perwindow/`: store individual folder for each window, including alignment and tree (if possible)
    - `wsize.perwindowsum`: summary table that shows if windows are informative or not according to `min_informative_sites` and `IQ-Tree2`
- `filtered/`: similar with `perwindow/` but with uninformative windows filtered out
- `step.aic.sum`: table that shows the delta AIC and chromosomal positions for each window(s) pair

### <a id="hmm">HMM</a>
*(to be implemented later)*

## <a id="datasets">Datasets</a>

### <a id="edelman">Genome of *erato-sara Heliconius* butterflies from <a href="https://doi.org/10.1126/science.aaw2090">Edelman et al. (2019)</a></a>
The dataset consists of 21 chromosomes of six *erato-sara Heliconius* species (*H. demeter, H. sara, H. telesiphe, H. hecalesia, H. himera, H. erato*) and one *H. melpomene* as outgroup.

#### Data Preparation
Following the original publication, we extract FASTA alignments from HAL by following these steps:
1. Convert HAL to MAF format using `hal2maf`
2. Retrieve single-copy MAF blocks using `getSingleCopy.py`
3. Sort the taxa of each MAF block using `maf-sort.sh`
4. Convert MAF blocks to FASTA alignments using `msa_view`
5. Concatenate the FASTA alignments using `msa_view`
6. (Additional step) Trim the leading and trailing positions of each chromosome based on the `gaps_threshold` value (see below).
7. (Additional step) Trim all positions of each chromosome based on the `gaps_threshold` value (see below).

| Parameters       | Definition                                                                                                      |
| ---------------- | --------------------------------------------------------------------------------------------------------------- |
| `outdir`         | Output directory                                                                                                |
| `thread`         | Number of threads for parallelization                                                                           |
| `redo`           | If `FALSE`, skip analysis if output files exist; if `TRUE`, overwrite previous results                          |
| `fn_hal`         | HAL file from Edelman et al. (2019), available <a href="https://doi.org/10.5061/dryad.b7bj832">here</a>         |
| `fn_refseq`      | List of reference sequence for HAL to MAF conversion, available <a href="edelman_etal_2019/refseq.txt">here</a> |
| `dir_hal2maf`    | Directory for `hal2maf` executable                                                                              |
| `dir_singleCopy` | Directory for `singleCopy.py` executable                                                                        |
| `dir_mafsort`    | Directory for `maf-sort.sh` executable                                                                          |
| `dir_msaview`    | Directory for `msa_view` executable                                                                             |
| `gaps_threshold` | Maximum gap proportion per site                                                                                 |

> **Important Notes** <br>
In `edelman_etal_2019/files/refseq.txt`, I comment out `Herato_chr17_3` as it does not contain single-copy MAF block, thus resulting in zero-length FASTA alignment.

#### Output
Running the code will create individual folder for each chromosome (i.e., chr1 to chr21). Each folder will contain the following subfolders:
- `raw/`: store the raw MAF blocks
- `singleCopy/`: store single-copy MAF blocks
- `sorted/`: store sorted MAF blocks
- `fasta/`: store FASTA alignment per MAF block
    - `concatenation/`: store raw and filtered concatenated FASTA alignments

```
outdir/
├── chr1/
│   ├── raw/
│   │   └── chr1_1.maf
│   ├── singleCopy/
│   │   └── chr1_1_singleCopy.maf
│   ├── sorted/
│   │   └── chr1_1_sorted.maf
│   └── fasta/
│       ├── chr1_1.fa
│       └── concatenation/
│           ├── chr1_concat.fa
│           ├── chr1_concat_nogaps.fa
│           └── chr1_concat_filtered.fa
...
├── chr21/
│   └── ...
├── edelman.html
└── edelman.log 
```

#### Time complexity
Running `run_all.R` using a server with `Intel(R) Xeon(R) CPU E5-2690 v4 @2.60GHz` and `Ubuntu 20.04.5 LTS`, the time required to run each step is as follows:
- Data preparation and filtering: ~1.6 hours with 50 threads
- Non-overlapping window analysis: up to 4.4 hours for the longest chromosome (`chr10`)
- Stepwise non-overlapping windows: `N/A`

#### Folder structure
<i>Non-overlapping window</i>
```
prefix/
├── chr1/
│   ├── chr1.fa
│   ├── chr1.log
│   ├── chr1.html
│   └── windows/
│       ├── 10000/
│       │   ├── alignment/
│       │   │   ├── window_01.fa
│       │   │   ├── window_02.fa
│       │   │   ├── window_03.fa
│       │   │   ├── window_04.fa
│       │   │   └── window_05.fa
│       │   ├── summary/
│       │   │   ├── prefix.atsum
│       │   │   ├── prefix.cnsum
│       │   │   └── prefix.topsum
│       │   └── trees/
│       │       ├── prefix.best_model.nex
│       │       ├── prefix.ckp.gz
│       │       ├── prefix.iqtree
│       │       ├── prefix.log
│       │       ├── prefix.parstree
│       │       └── prefix.treefile
│       ...
│       ├── 50000/
│       │   └── ...
│       ├── prefix.sum
│       ├── prefix.aicc.tiff
│       ├── prefix.aic.tiff
│       ├── prefix.bic.tiff
│       ├── prefix.topdist
│       └── prefix.topdist.tiff
...
└── chr21/
    └── ...
```

<i>Stepwise non-overlapping window</i>
```
prefix/
├── chr1/
│   ├── 1/
│   │   ├── 64000/
│   │   │   ├── perwindow/
│   │   │   │   ├── window_0001/
│   │   │   │   │   ├── window_01.fa
│   │   │   │   │   ├── window_01.fa.log
│   │   │   │   │   ├── window_01.fa.iqtree (if possible)
│   │   │   │   │   └── window_01.fa.treefile (if possible)
│   │   │   │   ...
│   │   │   │   ├── window_0262/
│   │   │   │   │   └── ...
│   │   │   │   └── 64000.perwindowsum
│   │   │   ├── filtered/
│   │   │   │   └── ...
│   │   │   └── 64000.uqtops
│   │   ├── 32000
│   │   │   └── ...
│   │   └── 1.aic.sum
│   ...
│   ├── 9/
│   │   └── ...
│   └── chr1.log
...
└── chr21/
    └── ...
```

---
*Last update: 06 November 2023 by Jeremias Ivan*