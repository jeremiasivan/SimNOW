# SimNOW

**SimNOW (Simulation Non-Overlapping Windows)** is an R pipeline to assess the accuracy of non-overlapping windows analysis in simulated phylogenetic studies. It consists of three main steps: sequence simulation, non-overlapping window analysis, and HMM.

## Table of Content
- <a href="#prereqs">Prerequisites</a>
- <a href="#genpipe">General Pipeline</a>
- <a href="#timecom">Time Complexity</a>
- <a href="#emps">Empirical Analyses</a>
- <a href="#refs">References</a>

## <a id="prereqs">Prerequisites</a>
This pipeline requires several software and R packages to run. All software have to be executable, while the R packages should be installed either in your local directory or virtual environment. We recommend you to use environment management system (e.g. `conda`) to install the packages, but you can also use `install.packages()` built-in function in R or RStudio.

### Software
- <a href="http://home.uchicago.edu/~rhudson1/source/mksamples.html">ms</a>
- <a href="http://www.iqtree.org/doc/AliSim">AliSim</a>
- <a href="http://www.iqtree.org">IQ-TREE 2</a>
- <a href="http://www.iqtree.org/doc/Complex-Models#multitree-models">MAST</a>-HMM
- <a href="https://bioinf.shenwei.me/seqkit/">Seqkit</a> (can also be installed via <a href="https://anaconda.org/bioconda/seqkit">Anaconda</a>)

### R packages
|    Name    |                               CRAN                               |                             Anaconda                             |
| ---------- |:----------------------------------------------------------------:|:----------------------------------------------------------------:|
| ape        | <a href="https://cran.r-project.org/package=ape">Link</a>        | <a href="https://anaconda.org/conda-forge/r-ape">Link</a>        |
| data.table | <a href="https://cran.r-project.org/package=data.table">Link</a> | <a href="https://anaconda.org/conda-forge/r-data.table">Link</a> |
| dplyr      | <a href="https://cran.r-project.org/package=dplyr">Link</a>      | <a href="https://anaconda.org/conda-forge/r-dplyr">Link</a>      |
| doSNOW     | <a href="https://cran.r-project.org/package=doSNOW">Link</a>     | <a href="https://anaconda.org/conda-forge/r-dosnow">Link</a>     |
| forcats    | <a href="https://cran.r-project.org/package=forcats">Link</a>    | <a href="https://anaconda.org/conda-forge/r-forcats">Link</a>    |
| ggplot2    | <a href="https://cran.r-project.org/package=ggplot2">Link</a>    | <a href="https://anaconda.org/conda-forge/r-ggplot2">Link</a>    |
| log4r      | <a href="https://cran.r-project.org/package=log4r">Link</a>      | <a href="https://anaconda.org/conda-forge/r-log4r">Link</a>      |
| rmarkdown  | <a href="https://cran.r-project.org/package=rmarkdown">Link</a>  | <a href="https://anaconda.org/conda-forge/r-rmarkdown">Link</a>  |
| seqinr     | <a href="https://cran.r-project.org/package=seqinr">Link</a>     | <a href="https://anaconda.org/conda-forge/r-seqinr">Link</a>     |
| stringr    | <a href="https://cran.r-project.org/package=stringr">Link</a>    | <a href="https://anaconda.org/conda-forge/r-stringr">Link</a>    |
| viridis    | <a href="https://cran.r-project.org/package=viridis">Link</a>    | <a href="https://anaconda.org/conda-forge/r-viridis">Link</a>    |

## <a id="genpipe">General Pipeline</a>
1. **Clone the Git repository** <br>
    ```
    git clone git@github.com:jeremiasivan/SimNOW.git
    ```

2. **Install the prerequisites** <br>
    > **Important Notes** <br>
    `AliSim` and `MAST-HMM` are integrated into `IQ-TREE 2`, but they are available in different releases. As a result, you need to individually download the latest version of `IQ-TREE 2` for both.

    Please download `ms` and `IQ-TREE 2` from the links above. For `Seqkit` and `R` packages, I prefer to download them from `Anaconda` as below.

    - Setting up conda environment with R
        ```
        conda create -n simnow -c conda-forge r-base
        conda activate simnow
        ```
    -  Installing R packages
        ```
        conda install seqkit
        conda install -c conda-forge r-ape r-data.table r-dplyr r-dosnow r-forcats r-ggplot2 r-log4r r-rmarkdown r-seqinr r-stringr r-viridis
        ```
        Notes: Here, I list all packages that are used in the pipeline, which is not ideal if you have a lot of packages. In that case, you can start by installing <a href="https://anaconda.org/conda-forge/r-essentials">`r-essentials`</a> which includes commonly-used packages in R.

3. **Update the parameters in the file** <br>
    Please refer to `codes/README.md` for the details of each parameter and which files to be updated. 

4. **Run the code file** <br>
    For running individual steps:
    ```
    Rscript -e "rmarkdown::render('~/SimNOW/codes/1_sequence_simulation/simulation.Rmd')"
    Rscript -e "rmarkdown::render('~/SimNOW/codes/2_non_overlapping_window/1_main.Rmd')"
    Rscript -e "rmarkdown::render('~/SimNOW/codes/3_hmm/1_main.Rmd')"
    Rscript -e "rmarkdown::render('~/SimNOW/codes/4_all_runs_summary/summary_all.Rmd')"
    ```

    For running the whole pipeline:
    ```
    Rscript ~/SimNOW/codes/run_all.R
    ```

    In UNIX-based operating systems (e.g., Linux and MacOS), it is advisable to use `nohup` or `tmux` to run the whole pipeline. For Windows, you can use `start`, but I have never tried it before. 

## <a id="timecom">Time Complexity</a>
Running the whole pipeline with the second scenario (see <a href="/codes/README.md#example">here</a>) using a server with `Intel(R) Xeon(R) CPU E5-2690 v4 @2.60GHz` and `Ubuntu 20.04.5 LTS`, the time required to run each step per replicate is as follows:
- Sequence simulation: up to 1 hour for `ms_r = 2000`
- Non-overlapping windows: ~2.25 hours
- HMM: up to 5 minutes for `ms_r = 2000`

## <a id="emps">Empirical Analyses</a>
In order to run non-overlapping windows on empirical datasets, please refer to <a href="codes/5_empirical_data/README.md">`codes/5_empirical_data/README.md`</a>. There are additional software and R packages required to run the pipeline, specific for each dataset:

### Edelman et al. (2019)
#### Software
- `hal2maf` from <a href="https://github.com/ComparativeGenomicsToolkit/hal">`HAL Toolkit`</a>
- `getSingleCopy.py` from <a href="https://doi.org/10.5281/zenodo.3401692">Edelman et al. (2019)</a>
- `maf-sort.sh` from <a href="https://github.com/UCSantaCruzComputationalGenomicsLab/last">`last`</a>
- `msa_view` from <a href="http://compgen.cshl.edu/phast/">`PHAST`</a>

#### R packages
|    Name    |                               CRAN                               |                             Anaconda                             |
| ---------- |:----------------------------------------------------------------:|:----------------------------------------------------------------:|
| numbers    | <a href="https://cran.r-project.org/package=numbers">Link</a>    | <a href="https://anaconda.org/conda-forge/r-numbers">Link</a>    |
| tidyr      | <a href="https://cran.r-project.org/package=tidyr">Link</a>      | <a href="https://anaconda.org/conda-forge/r-tidyr">Link</a>      |

---
## <a id="refs">References</a>
1. Hudson, R. (<a href="https://doi.org/10.1093/bioinformatics/18.2.337">2002</a>). **Generating samples under a Wright–Fisher neutral model of genetic variation**. *Bioinformatics*, *18*(2), 337–338.

2. Ly-Trong, N., et al. (<a href="https://doi.org/10.1093/molbev/msac092">2022</a>). **AliSim: A Fast and Versatile Phylogenetic Sequence Simulator for the Genomic Era**. *Molecular Biology and Evolution*, *39*(5), msac092.

3. Minh, B.Q., et al. (<a href="https://doi.org/10.1093/molbev/msaa015">2020</a>). **IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era**. *Molecular Biology and Evolution*, *37*(5), 1530–1534.
 
4. Wong, T., et al. (<a href="https://doi.org/10.1101/2022.10.06.511210">2022</a>). **MAST: Phylogenetic Inference with Mixtures Across Sites and Trees**. *bioRxiv*.

5. Shen, W., et al. (<a href="https://doi.org/10.1371/journal.pone.0163962">2016</a>). **SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation**. *PLoS ONE*, *11*(10), e0163962.

6. Edelman, et al. (<a href="https://doi.org/10.1126/science.aaw2090">2019</a>). **Genomic architecture and introgression shape a butterfly radiation**. *Science*, *366*(6465), 594–599.

7. Hickey, et al. (<a href="https://doi.org/10.1093/bioinformatics/btt128">2013</a>). **HAL: a hierarchical format for storing and analyzing multiple genome alignments**. *Bioinformatics*, *29*(10), 1341–1342.

8. Hubisz, et al. (<a href="https://doi.org/10.1093/bib/bbq072">2010</a>). **PHAST and RPHAST: phylogenetic analysis with space/time models**. *Briefings in Bioinformatics*, *12*(1), 41–51.

---
*Last update: 29 August 2023 by Jeremias Ivan*