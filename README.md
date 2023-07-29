# SimNOW

**SimNOW (Simulation Non-Overlapping Windows)** is an R pipeline to assess the accuracy of non-overlapping windows analysis in simulated phylogenetic studies.

---
## Prerequisites
This pipeline requires several software and R packages to run. 

### List of software:
- <a href="http://home.uchicago.edu/~rhudson1/source/mksamples.html">ms</a>
- <a href="http://www.iqtree.org/doc/AliSim">AliSim</a>
- <a href="http://www.iqtree.org">IQ-TREE 2</a>
- <a href="http://www.iqtree.org/doc/Complex-Models#multitree-models">MAST</a>

### List of R package:
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

---
## Notes
This pipeline consists of **sequence simulation**, **non-overlapping window analysis**, and **HMM**. Further details are available in `codes` folder.

---
## References
1. Hudson, R. (<a href="https://doi.org/10.1093/bioinformatics/18.2.337">2002</a>). **Generating samples under a Wright–Fisher neutral model of genetic variation**. *Bioinformatics*, *18*(2), 337–338.

2. Ly-Trong, N., et al. (<a href="https://doi.org/10.1093/molbev/msac092">2022</a>). **AliSim: A Fast and Versatile Phylogenetic Sequence Simulator for the Genomic Era**. *Molecular Biology and Evolution*, *39*(5), msac092.

3. Minh, B.Q., et al. (<a href="https://doi.org/10.1093/molbev/msaa015">2020</a>). **IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era**. *Molecular Biology and Evolution*, *37*(5), 1530–1534.
 
4. Wong, T., et al. (<a href="https://doi.org/10.1101/2022.10.06.511210">*in press*</a>). **MAST: Phylogenetic Inference with Mixtures Across Sites and Trees**. *bioRxiv*.

---
*Last update: 28 July 2023 by Jeremias Ivan*