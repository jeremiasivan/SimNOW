# SimNOW

**SimNOW (Simulation Non-Overlapping Windows)** is an R pipeline to assess the accuracy of non-overlapping windows approach in simulated phylogenetic studies. It starts with sequence simulation, followed by non-overlapping window analyses using fixed window sizes and/or variable window sizes. It is mainly developed and tested using MacOS and Linux, so there might be incompatibilities using Windows.

> [!NOTE]
> The stepwise non-overlapping window analyses have been moved to <a href="https://github.com/jeremiasivan/StepwiseNOW"><b>StepwiseNOW (Stepwise Non-Overlapping Windows)</b></a>. The codes are also archived at `archive-empirical` branch.

> [!NOTE]
> This repository stores the R codes for simulating and running non-overlapping windows on simulated chromosomes. If you have your own chromosome alignments, please use <a href="https://github.com/jeremiasivan/PhyloNOW"><b>PhyloNOW (Phylogenomic Non-Overlapping Windows)</b></a> to run non-overlapping windows with variable window sizes, or <a href="https://github.com/jeremiasivan/StepwiseNOW"><b>StepwiseNOW (Stepwise Non-Overlapping Windows)</b></a> to run non-overlapping windows with fixed window size (which has been shown to perform poorly compared to using variable window sizes).

**If you use SimNOW, please cite as:**
```
J. Ivan, P. Frandsen, R. Lanfear. (2026). Selecting a Window Size for Phylogenomic Analyses of Whole Genome Alignments Using AIC. Systematic Biology, 75(1), 100–114. doi:10.1093/sysbio/syaf053
```

## Table of Content
- <a href="#prereqs">Prerequisites</a>
- <a href="#genpipe">General Pipeline</a>
- <a href="#inout">Input and Output Files</a>
- <a href="#refs">References</a>

## <a id="prereqs">Prerequisites</a>
SimNOW requires several software and R packages to run. We recommend you to use environment management system (e.g. `conda`) to install the prerequisites, but you can also use `install.packages()` built-in function in R or RStudio.

### Software
|    Name    |                                  Website                                    |                             Anaconda                             |
| ---------- |:---------------------------------------------------------------------------:|:----------------------------------------------------------------:|
| IQ-TREE    | <a href="http://www.iqtree.org">Link</a>                                    | <a href="https://anaconda.org/bioconda/iqtree">Link</a>          |
| ms         | <a href="http://home.uchicago.edu/~rhudson1/source/mksamples.html">Link</a> | <a href="https://anaconda.org/bioconda/ms">Link</a>              |
| SeqKit     | <a href="https://bioinf.shenwei.me/seqkit/">Link</a>                        | <a href="https://anaconda.org/bioconda/seqkit">Link</a>          |


### R packages
|    Name    |                               CRAN                               |                             Anaconda                             |
| ---------- |:----------------------------------------------------------------:|:----------------------------------------------------------------:|
| ape        | <a href="https://cran.r-project.org/package=ape">Link</a>        | <a href="https://anaconda.org/conda-forge/r-ape">Link</a>        |
| data.table | <a href="https://cran.r-project.org/package=data.table">Link</a> | <a href="https://anaconda.org/conda-forge/r-data.table">Link</a> |
| doSNOW     | <a href="https://cran.r-project.org/package=doSNOW">Link</a>     | <a href="https://anaconda.org/conda-forge/r-dosnow">Link</a>     |
| kableExtra | <a href="https://cran.r-project.org/package=kableExtra">Link</a> | <a href="https://anaconda.org/conda-forge/r-kableextra">Link</a> |
| log4r      | <a href="https://cran.r-project.org/package=log4r">Link</a>      | <a href="https://anaconda.org/conda-forge/r-log4r">Link</a>      |
| optparse   | <a href="https://cran.r-project.org/package=optparse">Link</a>   | <a href="https://anaconda.org/conda-forge/r-optparse">Link</a>   |
| rmarkdown  | <a href="https://cran.r-project.org/package=rmarkdown">Link</a>  | <a href="https://anaconda.org/conda-forge/r-rmarkdown">Link</a>  |
| seqinr     | <a href="https://cran.r-project.org/package=seqinr">Link</a>     | <a href="https://anaconda.org/conda-forge/r-seqinr">Link</a>     |
| tidyverse  | <a href="https://cran.r-project.org/package=tidyverse">Link</a>  | <a href="https://anaconda.org/conda-forge/r-tidyverse">Link</a>  |
| viridis    | <a href="https://cran.r-project.org/package=viridis">Link</a>    | <a href="https://anaconda.org/conda-forge/r-viridis">Link</a>    |
| yaml       | <a href="https://cran.r-project.org/package=yaml">Link</a>       | <a href="https://anaconda.org/conda-forge/r-yaml">Link</a>       |

## <a id="genpipe">General Pipeline</a>
1. **Clone the Git repository** <br>
    ```
    git clone git@github.com:jeremiasivan/SimNOW.git
    ```

2. **Install the prerequisites** <br>
    - Create a new conda environment
        ```
        conda create -n simnow
        conda activate simnow
        ```
    -  Installing prerequisites
        ```
        conda install -c conda-forge r-ape r-data.table r-doSNOW r-kableExtra r-log4r r-optparse r-rmarkdown r-seqinr r-tidyverse r-viridis r-yaml bioconda::iqtree bioconda::ms bioconda::seqkit
        ```

3. **Update the parameters in `config.yaml`** <br>

4. **Run SimNOW** <br>
    ```
    Rscript run_pipeline.R --config config.yaml
    Rscript run_pipeline.R --config config.yaml --redo
    ```

    In UNIX-based operating systems (e.g., Linux and MacOS), it is advisable to use `nohup` or `tmux` to run the whole pipeline. For Windows, you can use `psmux`.

## <a id="inout">Input and Output Files</a>
Please see <a href="./codes/README.md">README.md</a> for more details about input and output files required for running SimNOW, including examples from publications.

---
## <a id="refs">References</a>
1. Ly-Trong, N., et al. (<a href="https://doi.org/10.1093/molbev/msac092">2022</a>). **AliSim: A Fast and Versatile Phylogenetic Sequence Simulator for the Genomic Era**. *Molecular Biology and Evolution*, *39*(5), msac092.

2. Minh, B.Q., et al. (<a href="https://doi.org/10.1093/molbev/msaa015">2020</a>). **IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era**. *Molecular Biology and Evolution*, *37*(5), 1530–1534.

3. Hudson, R. (<a href="https://doi.org/10.1093/bioinformatics/18.2.337">2002</a>). **Generating samples under a Wright–Fisher neutral model of genetic variation**. *Bioinformatics*, *18*(2), 337–338.

4. Shen et al. (<a href="https://doi.org/10.1371/journal.pone.0163962">2016</a>). **SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation**. *PLOS ONE*, *11*(10), e0163962.

5. Paradis, E., et al. (<a href="https://doi.org/10.1093/bioinformatics/btg412">2004</a>). **APE: Analyses of Phylogenetics and Evolution in R language**. *Bioinformatics*, *20*(2), 289-290.

6. Barrett, T., et al. (<a href="https://doi.org/10.32614/CRAN.package.data.table">2026</a>). **data.table: Extension of 'data.frame'**. *R package*.

7. Daniel, F. (<a href="https://cran.r-project.org/package=doSNOW">2022</a>). **doSNOW: Foreach Parallel Adaptor for the 'snow' Package**. *R package*.

8. Zhu, H., et al. (<a href="https://cran.r-project.org/package=kableExtra">2024</a>). **kableExtra: Construct Complex Table with 'kable' and Pipe Syntax**. *R package*.

9. White, J.M., & Jacobs, A. (<a href="https://doi.org/10.32614/CRAN.package.log4r">2024</a>). **log4r: A Fast and Lightweight Logging System for R, Based on 'log4j'**. *R package*.

10. Davis, T.L. (<a href="https://doi.org/10.32614/CRAN.package.optparse">2026</a>). **optparse: Command Line Option Parser**. *R package*.

11. Allaire, J.J., et al. (<a href="https://doi.org/10.32614/CRAN.package.rmarkdown">2026</a>). **rmarkdown: Dynamic Documents for R**. *R package*.

12. Charif, D., et al. (<a href="https://cran.r-project.org/package=seqinr">2024</a>). **seqinr: Biological Sequences Retrieval and Analysis**. *R package*.

13. Wickham, H., et al. (<a href="https://doi.org/10.21105/joss.01686">2019</a>). **Welcome to the tidyverse**. *Journal of Open Source Software*, *4*(43), 1686.

14. Garnier, S., et al. (<a href="https://cran.r-project.org/package=viridis">2024</a>). **viridis: Colorblind-Friendly Color Maps for R**. *R package*.

15. Stephens, J., et al. (<a href="https://doi.org/10.32614/CRAN.package.yaml">2025</a>). **yaml: Methods to Convert R Data to YAML and Back**. *R package*.

16. Anthropic. (<a href="https://claude.ai/">2026</a>). Claude 4.6 Sonnet was used to generate `config.yaml` and `run_pipeline.R`. 

---
*Last update: 27 June 2026 by Jeremias Ivan*