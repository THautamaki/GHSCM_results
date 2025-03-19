
# Results of the article: Gene co-expression network estimation using GEM with GHS prior

This repository contains all codes and results files used to produce the
results for the article *Gene co-expression network estimation using
generalized expectation-maximization algorithm with graphical horseshoe
prior*.

## Original MATLAB codes

The original MATLAB codes of the GHS Gibbs sampler (GHS MCMC), GHS local
linear approximation (GHS LLA), GHS-like expectation conditional
maximization (GHS-like ECM or HSL ECM), and GHS-like MCMC (or HSL MCMC)
methods can be found in the following GitHub pages:

- <https://github.com/liyf1988/GHS>
  - Original author: Yunfan Li
- <https://github.com/sagarknk/GHS-LLA-codes>
  - Original author: Ksheera Sagar K. N., Purdue University
  - Codes modified and redistributed with his permission.
- <https://github.com/sagarknk/Graphical_HSL>
  - Original author: Ksheera Sagar K. N., Purdue University
  - Codes modified and redistributed with his permission.

## Real-world datasets

The CEU dataset is available on Sanger Institute’s website
<https://ftp.sanger.ac.uk/pub/genevar> and the DLBC dataset is available
on the website of the PRECISE framework
<https://mjha.shinyapps.io/PRECISE>.

## Structure of the repository

The files in the repository are organized in the main directories as
follows:

- Data
  - This directory is not stored in the GitHub as the simulated datasets
    can be generated at any time and real-world datasets are available
    publicly on the internet. The R code creates this when simulated
    datasets are generated. Please store the real-world datasets in the
    directories named “CEU_dataset” and “DLBC_dataset”.
- Figures
  - Contains all figures generated for the main article and
    supplementary material in own directories.
- MATLAB_files
  - The root contains some helper functions.
  - Each method have own directories where associated codes are.
  - Both real-world datasets also have own directories where are all
    MATLAB codes needed to analyse the datasets.
- R_files
  - Contains all codes used in the R.
- Results_files
  - The root contains the total times of the MATLAB methods.
  - The simulations with both problem dimensions $p = \{100, 200\}$ have
    own directories, where results are stored per method per directory.
  - Both real-world datasets have own directories containing all results
    files.

## Performing the simulation analyses

> [!NOTE\
> All provided scripts assume the working directory has the
> same directory structure as this repository. The scripts use relative
> paths to read and save files.

To run analyses using the R methods, we created the R script named
`Run_R_methods_and_combine_results.R`, the main file for the simulation
analysis. It is well commented so we do not give details here. It
installs all needed R packages, creates simulated datasets, runs R
methods, combines all results files, and finally prints results.

The simulated datasets used in the article can always be created using
the function `generate_datasets()` with no arguments (the default
arguments).

Each MATLAB method has `<method_name>_analysis.m` file to run
simulations. They are all well commented and can be run using MATLAB’s
`Run all sections` shortcut. Please note that the helper functions
should be placed into the `%userprofile%\Documents\MATLAB` folder (if a
Windows machine) so that MATLAB can find these.

Please note that if a similar system is used, as in the article, all
simulations take about 3.5 days to run.

## Performing the real-world dataset analyses

### The CEU dataset

We created R and MATLAB scripts to perform analyses of the CEU dataset.
We suggest creating transformed datasets first by running lines 7–32
from the R script `Analysis_of_CEU_dataset.R`. Then MATLAB methods can
be run from the script `Run_CEU_dataset_MATLAB_methods.m` MATLAB’s
`Run all sections` shortcut. Finally, rest of the lines from the R
script can be run. It performs the GHS GEM method, calculates numbers
used in the article and plot the estimated networks.

### The DLBC dataset

We created R and MATLAB scripts to perform analyses of the DLBC dataset.
We suggest creating transformed datasets first by running lines 8–35
from the R script `Analysis_of_DLBC_dataset.R`. Then MATLAB methods can
be run from the script `Run_DLBC_dataset_MATLAB_methods.m` MATLAB’s
`Run all sections` shortcut. Finally, the rest of the lines from the R
script can be run. It performs the GHS GEM method, calculates the
numbers used in the article, and plots the estimated networks.

## Supplementary materials

### Analysis of the kappa threshold

Figures S1, S2, and S3 in the supplementary material can be created
using the R script `Analysis_of_kappa_threshold.R`.

### Analysis of the initial values of the GHS GEM algorithm

Analysis of the initial values of the GHS GEM algorithm in the
supplementary material can be run using the R script
`Analysis_of_initial_values.R`.

### Parameter c in the global scale parameter formula

The R script `Analysis_of_c_parameter.R`.

### Results with the fastGHS method

We did not include the fastGHS method results in the main article, but
they can be reproduced using the same R script as used for the main
results (`Run_R_methods_and_combine_results.R`, lines from 102 forward).
We also did some additional analysis, which can be reproduced using the
R script `Additional_fastGHS_analysis.R`.
