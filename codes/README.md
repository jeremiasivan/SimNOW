# SimNOW

## General Pipeline
1. **Clone the Git repository**
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

## Analyses
> **Important Notes** <br>
In this pipeline, the results of one step might be used in the other steps. Please be careful about the *consistency* of the parameters and folder structure.

### Sequence simulation
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

**Example from the paper** <br>
In the paper, we want to simulate a *choromosome*. There are six parameters of interest from the list:
```
ms_params: "7 1 -T -I 7 1 1 1 1 1 1 1 -ej 2.13 2 1 -ej 5.60 3 1 -ej 7.07 4 1 -ej 6.67 5 6 -ej 8.13 6 1 -ej 15.47 7 1 -es 0.21 2 0.966 -ej 0.21 8 1 -es 0.21 1 0.796 -ej 0.21 9 2 -es 5.52 3 0.115 -ej 5.52 10 4 -es 6.84 6 0.778 -ej 6.84 11 4 -es 6.84 4 0.853 -ej 6.84 12 6"

ms_r: 0
ms_l: 1000000

alisim_model: "JC"
alisim_scale: 0.005
  
copy_gaps: FALSE
```

- `ms_params`: the first two parameters refer to the number of taxa (`7`) and chromosome (`1`). `-T` sets output as trees. `-I` specifies 7 populations with 1 sample from each population, which gives us the *multispecies coalescent model*. The next `-ej` flags specify the speciation events with three parameters: time in `4Ne` generation, and the two participating taxa. The `-es` flag reflects introgression event with three parameters: time in `4Ne generation`, recipient taxon, and inheritance probability (which is 1-introgression rate), followed by `-ej` with a newly split population and the donor taxon, respectively.
- `ms_r`: recombination rate in `ms`. Higher rate will result in more frequent tree switching across the chromosome.
- `ms_l`: total length of the chromosome.
- `alisim_model`: we used `JC` to simplify the simulation.
- `alisim_scale`: as `ms` uses `4Ne` generation as time unit, while `AliSim` uses substitution-per-site, we have to rescale the branch lengths by multiplying them with a scaling factor. This will affect the number of informative sites of the alignment.
- `copy_gaps`: we set this as `FALSE` to simplify the simulation.

**Output** <br>
Running the code will create a new folder called `simulation` with the following files:
- `prefix.ms`: output `ms` file
- `prefix.mstrees`: table of locus boundaries and trees
- `prefix.mstops`: table of locus boundaries and topologies
- `prefix.uqtops`: list of unique topology
- `prefix.nex`: partition file based on `ms` output trees
- `prefix.multitree`: list of locus trees from `prefix.mstrees`
- `prefix.fa`: FASTA alignment from `AliSim`
- `prefix_nogaps.fa`: FASTA alignment without gaps (only when `only_gaps` == `TRUE`)

### Non-overlapping window analysis
In this step, we use the simulated alignment to perform non-overlapping window analysis and generate the summary statistics.

### HMM
In this step, we use MAST-HMM to refine the breakpoints between windows. We simplify the process by running MAST-HMM on two consecutive windows if the windows have different topology.

## Notes
If you want to run more than one analysis at the same time, you can run `run_all.R`.

---
*Last update: 27 July 2023 by Jeremias Ivan*