# SimNOW

**SimNOW (Simulation Non-Overlapping Windows)** is an R pipeline to assess the accuracy of non-overlapping windows analysis in simulated phylogenetic studies.

Software incorporated in the pipeline:
- <a href="http://home.uchicago.edu/~rhudson1/source/mksamples.html">ms</a>
- <a href="http://www.iqtree.org/doc/AliSim">AliSim</a>
- <a href="http://www.iqtree.org">IQ-TREE 2</a>
- <a href="http://www.iqtree.org/doc/Complex-Models#multitree-models">MAST</a>+HMM

---
### Notes
1. Currently, you can only set `ms` recombination rate to vary between simulations via `wrapper/nw_main.R`
2. Multiple MAST models should be manually specified in `wrapper/mh_main.R` (complex models are unsupported)
3. Treelikeness test has *not* been tested

---
### Considerations
- Allow users to choose which analysis to run
- Automate the selection of MAST model based on the number of topologies
- Fix the summary report for the whole pipeline, including HTML template
- Parallelization, particularly for `nowsummary` block in `rmd/nw_run.Rmd`

---
### Future Work(s)
- Incorporate other methods to find correlates with non-overlapping window accuracy (e.g., information criteria, <a href="https://doi.org/10.1186/s12862-016-0837-3">phylogenetic information</a>, <a href="https://doi.org/10.1093/molbev/msaa106">gene</a> and <a href="https://doi.org/10.1093/bioinformatics/btac741">site concordance factor</a>, <a href="https://doi.org/10.1101/2021.02.16.431544">treelikeness</a>)

---
### References
1. Hudson, R. (<a href="https://doi.org/10.1093/bioinformatics/18.2.337">2002</a>). **Generating samples under a Wright–Fisher neutral model of genetic variation**. *Bioinformatics*, *18*(2), 337–338.

2. Ly-Trong, N., et al. (<a href="https://doi.org/10.1093/molbev/msac092">2022</a>). **AliSim: A Fast and Versatile Phylogenetic Sequence Simulator for the Genomic Era**. *Molecular Biology and Evolution*, *39*(5), msac092.

3. Minh, B.Q., et al. (<a href="https://doi.org/10.1093/molbev/msaa015">2020</a>). **IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era**. *Molecular Biology and Evolution*, *37*(5), 1530–1534.
 
4. Wong, T., et al. (<a href="https://doi.org/10.1101/2022.10.06.511210">*in press*</a>). **MAST: Phylogenetic Inference with Mixtures Across Sites and Trees**. *bioRxiv*.