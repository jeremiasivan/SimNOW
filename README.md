# SimNOW

**SimNOW (Simulation Non-Overlapping Windows)** is an R pipeline to assess the accuracy of non-overlapping windows analysis in simulated phylogenetic studies.

Software incorporated in the pipeline:
- <a href="http://home.uchicago.edu/~rhudson1/source/mksamples.html">ms</a>
- <a href="http://www.iqtree.org/doc/AliSim">AliSim</a>
- <a href="http://www.iqtree.org">IQ-TREE 2</a>
- <a href="http://www.iqtree.org/doc/Complex-Models#multitree-models">MAST</a>

---
### Notes
1. Currently, you can only set `ms` recombination rate to vary between simulations via `wrapper/nw_main.R`
2. Multiple MAST models should be manually specified in `wrapper/mh_main.R`
3. Some parts of the codes have not been fully tested

---
### Considerations
- Fix naming convention for file extension
- Allow users to choose which analysis to run
- Find a better way to _redo_ the analysis (e.g., delete all files within the directory, create a new directory, etc.)
- Automate the selection of MAST model based on the number of topologies
- Run MASTHMM based on unique set of input topologies instead of windows
- Fix the summary report for the whole pipeline, including HTML template

---
### Future Work(s)
- Incorporate other methods to assess (and improve) non-overlapping window accuracy

---
### References
1. Hudson, R. (<a href="https://doi.org/10.1093/bioinformatics/18.2.337">2002</a>). **Generating samples under a Wright–Fisher neutral model of genetic variation**. *Bioinformatics*, *18*(2), 337–338.

2. Ly-Trong, N., et al. (<a href="https://doi.org/10.1093/molbev/msac092">2022</a>). **AliSim: A Fast and Versatile Phylogenetic Sequence Simulator for the Genomic Era**. *Molecular Biology and Evolution*, *39*(5), msac092.

3. Minh, B.Q., et al. (<a href="https://doi.org/10.1093/molbev/msaa015">2020</a>). **IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era**. *Molecular Biology and Evolution*, *37*(5), 1530–1534.
 
4. Wong, T., et al. (<a href="https://doi.org/10.1101/2022.10.06.511210">*in press*</a>). **MAST: Phylogenetic Inference with Mixtures Across Sites and Trees**. *bioRxiv*.