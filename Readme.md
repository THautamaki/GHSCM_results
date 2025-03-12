
# Results of the article: Gene co-expression network estimation using GEM with GHS prior

This repository contains all codes used to produce the results for the
article *Gene co-expression network estimation using generalized
expectation-maximization algorithm with graphical horseshoe prior*.

## Structure of the repository

The files in the repository are organized in the directories as follows:

- R_files
  - Contains all codes used in the R.
- MATLAB_files
  - The root contains some helper functions.
    - CEU_dataset
      - Contains all MATLAB codes needed to analyze the CEU dataset.
    - DLBC_dataset
      - Contains all MATLAB codes needed to analyze the DLBC dataset.
    - GHS_LLA
      - Contains all MATLAB codes needed to run simulations using the
        GHS local linear approximation (LLA) algorithm with Laplace
        representation.
    - GHS_MCMC
      - Contains all MATLAB codes needed to run simulations using the
        Gibbs sampler with the GHS prior.
    - HSL_ECM
      - Contains all MATLAB codes needed to run simulations using the
        expectation conditional maximization (ECM) algorithm with the
        GHS-like prior.
    - HSL_MCMC
      - Contains all MATLAB codes needed to run simulations using the
        MCMC algorithm with the GHS-like prior.
- Results_files
  - The root contains the total times of the MATLAB methods.
  - CEU_dataset
  - DLBC_dataset
  - Contains directories for the simulations with both problem
    dimensions $p = \{100, 200\}$ and real-world datasets containing all
    results files. The root contains the total times of the MATLAB
    methods.
