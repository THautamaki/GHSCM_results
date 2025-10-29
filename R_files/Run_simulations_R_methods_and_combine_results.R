# Install all required packages.
if(!require("BDgraph", quietly = TRUE)) {
  install.packages("BDgraph")
}
if(!require("huge", quietly = TRUE)) {
  install.packages("huge")
}
if(!require("GHSCM", quietly = TRUE)) {
  devtools::install_github("THautamaki/GHSCM")
}
if(!require("beam", quietly = TRUE)) {
  devtools::install_version("beam", version = "2.0.2")
}
if(!require("fastGHS", quietly = TRUE)) {
  devtools::install_github("Camiling/fastGHS")
}
if(!require("doParallel", quietly = TRUE)) {
  install.packages("doParallel")
}
if(!require("pulsar", quietly = TRUE)) {
  install.packages("pulsar")
}
if(!require("igraph", quietly = TRUE)) {
  install.packages("igraph")
}
if(!require("R.matlab", quietly = TRUE)) {
  install.packages("R.matlab")
}

# Load packages.
library(GHSCM)
library(beam)
library(huge)
library(fastGHS)
library(doParallel)
library(pulsar)

# Load functions needed for the analyses.
source("R_files/Generate_datasets.R")
source("R_files/Functions_for_R_methods_and_combining_results.R")

# Change this FALSE if data is not yet generated.
data_generated <- TRUE
if (!data_generated) {
  generate_datasets()
}

# Set all needed parameters.
sample_sizes <- c(120)
variable_numbers <- c(100, 200)
structures <- c("random", "bdgraph_sf", "huge_sf", "hubs")
R_methods <- c("GHSCM", "GLASSO", "fastGHS")
MATLAB_methods <- c("GHS_MCMC", "GHS_LLA", "HSL_MCMC", "HSL_ECM")

# Run all R methods if not yet run.
simulations_done <- TRUE
if (!simulations_done) {
  all_results <- run_multiple_methods(R_methods, structures, sample_sizes, variable_numbers,
                                      save_results = TRUE)
}

# Combine results into one object.
# If results using R methods exist and saved, create empty list first.
all_results <- list()
# Read results from the results files.
all_results <- add_R_results(all_results, R_methods, structures, sample_sizes, variable_numbers)
all_results <- add_MATLAB_results(all_results, MATLAB_methods,
                                  structures, sample_sizes, variable_numbers)

# Define the methods and which order they will be printed. GHS CM has two results for random network
# structure.
random_methods <- c("GHSCM_p/2", "GHSCM", "GHS_MCMC", "fastGHS", "GHS_LLA", "HSL_MCMC", "HSL_ECM", "GLASSO")
other_methods <- c("GHSCM", "GHS_MCMC", "fastGHS", "GHS_LLA", "HSL_MCMC", "HSL_ECM", "GLASSO")

# Define scores which results will be printed.
scores <- c("MCC", "TPR", "FPR", "FDR", "f_norm_rel", "sl_omega", "time", "total_time")

# Print results.
print_results(scores, all_results, random_methods, "random", 120, 100)
print_results(scores, all_results, other_methods, "bdgraph_sf", 120, 100)
print_results(scores, all_results, other_methods, "huge_sf", 120, 100)
print_results(scores, all_results, other_methods, "hubs", 120, 100)

print_results(scores, all_results, random_methods, "random", 120, 200)
print_results(scores, all_results, other_methods, "bdgraph_sf", 120, 200)
print_results(scores, all_results, other_methods, "huge_sf", 120, 200)
print_results(scores, all_results, other_methods, "hubs", 120, 200)

### Create bodies of LaTeX tabels in the main paper.

# Define network estimation scores for LaTeX table in main paper (Table 2).
scores <- c("MCC", "TPR", "FPR", "FDR")

# Print bodies of LaTeX tabel.
create_latex_table(scores, all_results, random_methods, "random", 120, 100)
create_latex_table(scores, all_results, other_methods, "bdgraph_sf", 120, 100)
create_latex_table(scores, all_results, other_methods, "huge_sf", 120, 100)
create_latex_table(scores, all_results, other_methods, "hubs", 120, 100)

create_latex_table(scores, all_results, random_methods, "random", 120, 200)
create_latex_table(scores, all_results, other_methods, "bdgraph_sf", 120, 200)
create_latex_table(scores, all_results, other_methods, "huge_sf", 120, 200)
create_latex_table(scores, all_results, other_methods, "hubs", 120, 200)

# Define precision matrix estimation scores for LaTeX table in main paper (Table 3).
scores <- c("f_norm_rel", "sl_omega")

# Print bodies of LaTeX tabel.
create_latex_table(scores, all_results, random_methods, "random", 120, 100)
create_latex_table(scores, all_results, other_methods, "bdgraph_sf", 120, 100)
create_latex_table(scores, all_results, other_methods, "huge_sf", 120, 100)
create_latex_table(scores, all_results, other_methods, "hubs", 120, 100)

create_latex_table(scores, all_results, random_methods, "random", 120, 200)
create_latex_table(scores, all_results, other_methods, "bdgraph_sf", 120, 200)
create_latex_table(scores, all_results, other_methods, "huge_sf", 120, 200)
create_latex_table(scores, all_results, other_methods, "hubs", 120, 200)

# Create body of the run time table (Table 4).
create_latex_table("time", all_results, random_methods, "random", 120, 100)
create_latex_table("time", all_results, random_methods, "random", 120, 200)
create_latex_table("time", all_results, other_methods, "bdgraph_sf", 120, 100)
create_latex_table("time", all_results, other_methods, "bdgraph_sf", 120, 200)
create_latex_table("time", all_results, other_methods, "huge_sf", 120, 100)
create_latex_table("time", all_results, other_methods, "huge_sf", 120, 200)
create_latex_table("time", all_results, other_methods, "hubs", 120, 100)
create_latex_table("time", all_results, other_methods, "hubs", 120, 200)

create_latex_table("total_time", all_results, random_methods, "random", 120, 100)
create_latex_table("total_time", all_results, random_methods, "random", 120, 200)
create_latex_table("total_time", all_results, other_methods, "bdgraph_sf", 120, 100)
create_latex_table("total_time", all_results, other_methods, "bdgraph_sf", 120, 200)
create_latex_table("total_time", all_results, other_methods, "huge_sf", 120, 100)
create_latex_table("total_time", all_results, other_methods, "huge_sf", 120, 200)
create_latex_table("total_time", all_results, other_methods, "hubs", 120, 100)
create_latex_table("total_time", all_results, other_methods, "hubs", 120, 200)

#######
# Next lines are all for Appendix.

####
# GHS LLA runtimes and proportion of the tau's tuning times for Appendix D (Table D.3).
for (n in sample_sizes) {
  for (p in variable_numbers) {
    cat("n: ", n, ", p: ", p, "\n", sep = "")
    for (structure in structures) {
      tau_prop <- mean(all_results$GHS_LLA[[structure]][[paste0("n", n, "_p", p)]]$results$tau_prop)
      tau_time <- mean(all_results$GHS_LLA[[structure]][[paste0("n", n, "_p", p)]]$results$tau_time)
      alg_time <- mean(all_results$GHS_LLA[[structure]][[paste0("n", n, "_p", p)]]$results$alg_time)
      cat(format(structure, width = 10), " & ", round(tau_time), " & ", round(alg_time), " & ", round(tau_prop, 3)*100, "\n", sep = "")
    }
  }
}

####
# Create LaTeX table body of mean number of false positives for Appendix D (Table D.2).
create_false_positives_table(all_results, random_methods, structures, sample_sizes, variable_numbers)

####
# Create LaTex table body of mean number of sparsities for Appendix D (Table D.1).
create_sparsity_table(all_results, random_methods, structures, sample_sizes, variable_numbers)

# True sparsities
# p = 100
# BDgraph random
round(100 - (50 / (100 * (100 - 1) / 2) * 100), 2)
# Scale-free networks (both BDgraph and huge are the same)
round(100 - (99 / (100 * (100 - 1) / 2) * 100), 2)
# Hub
round(100 - (95 / (100 * (100 - 1) / 2) * 100), 2)

# p = 200
# BDgraph random
round(100 - (100 / (200 * (200 - 1) / 2) * 100), 2)
# Scale-free networks (both BDgraph and huge are the same)
round(100 - (199 / (200 * (200 - 1) / 2) * 100), 2)
# Hub
round(100 - (190 / (200 * (200 - 1) / 2) * 100), 2)
