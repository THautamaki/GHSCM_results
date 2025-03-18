
# Results of the article: Gene co-expression network estimation using GEM with GHS prior

This repository contains all codes used to produce the results for the
article *Gene co-expression network estimation using generalized
expectation-maximization algorithm with graphical horseshoe prior*.

## Original MATLAB codes

The original MATLAB codes for the GHS Gibbs sampler (GHS MCMC), GHS
local linear approximation (GHS LLA), GHS-like expectation conditional
maximization (GHS-like ECM or HSL ECM), and GHS-like MCMC (or HSL MCMC)
methods can be found in the following GitHub pages:

- <https://github.com/liyf1988/GHS>
- <https://github.com/sagarknk/GHS-LLA-codes>
- <https://github.com/sagarknk/Graphical_HSL>

## Real-world datasets

The CEU dataset is available on Sanger Institute’s website
<https://ftp.sanger.ac.uk/pub/genevar> and the DLBC dataset is available
on the website of the PRECISE framework
<https://mjha.shinyapps.io/PRECISE>.

## Structure of the repository

The files in the repository are organized in the main directories as
follows:

- Data
  - This directory is not stored in the GitHub as simulated datasets can
    be generated at any time and real-world datasets are available
    publicly on the internet. The R code creates this when simulated
    datasets are generated.
- Figures
  - Contains all figures generated for the main article and
    supplementary material in own directories.
- MATLAB_files
  - The root contains some helper functions.
  - Each method have own directories where associated codes are.
  - Both real-world datasets have also own directories where all MATLAB
    codes needed to analyse the datasets are.
- R_files
  - Contains all codes used in the R.
- Results_files
  - The root contains the total times of the MATLAB methods.
  - The simulations with both problem dimensions $p = \{100, 200\}$ have
    own directories, where results are stored per method per directory.
  - Both real-world datasets have own directories containing all results
    files.

## Performing the simulations

### Simulation analyses

To run analyses using R methods, we created the R script named
Run_R_methods_and_combine_results.R, which is main file for the
simulation analysis. It is well commented so we do not give details
here.

Each MATLAB method have <method_name>\_analysis.m file to run
simulations. They are all well commented, and they can be run using
MATLAB’s Run all sections shortcut. Please note that helper functions
should be placed into %userprofile%\Documents\MATLAB folder (if a
Windows machine) that MATLAB can found these.
