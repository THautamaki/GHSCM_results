# Install all required packages.
if(!require("BDgraph", quietly = TRUE)) {
  install.packages("BDgraph")
}
if(!require("huge", quietly = TRUE)) {
  install.packages("huge")
}
if(!require("GHSGEM", quietly = TRUE)) {
  devtools::install_github("THautamaki/GHSGEM")
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

# Load packages.
library(GHSGEM)
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
R_methods <- c("GHSGEM", "GLASSO", "beam")
MATLAB_methods <- c("GHS_MCMC", "GHS_LLA", "HSL_MCMC", "HSL_ECM")

# Run all R methods if not yet run.
all_results <- run_multiple_methods(R_methods, structures, sample_sizes, variable_numbers,
                                    save_results = TRUE)

# Combine results into one object.
# If results using R methods exist and saved, create empty list first.
all_results <- list()
# Read results from the results files.
all_results <- add_R_results(all_results, R_methods, structures, sample_sizes, variable_numbers)
all_results <- add_MATLAB_results(all_results, MATLAB_methods,
                                  structures, sample_sizes, variable_numbers)

# Define the methods and which order they will be printed. GHS GEM has two results for random network
# structure.
random_methods <- c("GHSGEM_p/2", "GHSGEM", "GHS_MCMC", "GHS_LLA", "HSL_MCMC", "HSL_ECM", "GLASSO", "beam")
other_methods <- c("GHSGEM", "GHS_MCMC", "GHS_LLA", "HSL_MCMC", "HSL_ECM", "GLASSO", "beam")

# Define scores which results will be printed.
scores <- c("MCC", "TPR", "FPR", "FDR", "f_norm_rel", "sl_omega", "time")

# Print results.
print_results(scores, all_results, random_methods, "random", 120, 100)
print_results(scores, all_results, other_methods, "bdgraph_sf", 120, 100)
print_results(scores, all_results, other_methods, "huge_sf", 120, 100)
print_results(scores, all_results, other_methods, "hubs", 120, 100)

print_results(scores, all_results, random_methods, "random", 120, 200)
print_results(scores, all_results, other_methods, "bdgraph_sf", 120, 200)
print_results(scores, all_results, other_methods, "huge_sf", 120, 200)
print_results(scores, all_results, other_methods, "hubs", 120, 200)

# Define scores for LaTeX tables in main paper.
scores <- c("MCC", "TPR", "FDR", "sl_omega")

# Print LaTeX tabels.
# create_latex_table(scores, all_results, random_methods, "random", 120, 100)
# create_latex_table(scores, all_results, other_methods, "bdgraph_sf", 120, 100)
# create_latex_table(scores, all_results, other_methods, "huge_sf", 120, 100)
# create_latex_table(scores, all_results, other_methods, "hubs", 120, 100)

# Table 1.
create_latex_table(scores, all_results, random_methods, "random", 120, 200)
create_latex_table(scores, all_results, other_methods, "bdgraph_sf", 120, 200)
create_latex_table(scores, all_results, other_methods, "huge_sf", 120, 200)
create_latex_table(scores, all_results, other_methods, "hubs", 120, 200)

#######
# Next lines are all for Supplementary Materials.

######
# Create full LaTeX tables for supplementary material.
# Define scores.
scores <- c("MCC", "TPR", "FPR", "FDR", "f_norm_rel", "sl_omega", "time")

# Print LaTeX tabels. First, p = 100.
# Supplementary Table S4.
create_latex_table(scores, all_results, random_methods, "random", 120, 100)
create_latex_table(scores, all_results, other_methods, "bdgraph_sf", 120, 100)
create_latex_table(scores, all_results, other_methods, "huge_sf", 120, 100)
create_latex_table(scores, all_results, other_methods, "hubs", 120, 100)

# Next, p = 200
# Supplementary Table S5.
create_latex_table(scores, all_results, random_methods, "random", 120, 200)
create_latex_table(scores, all_results, other_methods, "bdgraph_sf", 120, 200)
create_latex_table(scores, all_results, other_methods, "huge_sf", 120, 200)
create_latex_table(scores, all_results, other_methods, "hubs", 120, 200)

#######
# GHS LLA runtimes and proportion of the tau's tuning times.
# Supplementary Table S3.
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

######
# Create LaTeX table of mean number of false positives for supplementary material.
# Supplementary Table S6.
create_false_positives_table(all_results, random_methods, structures, sample_sizes, variable_numbers)


# Results of the fastGHS method with p = 100.
# Initialize paramaters.
sample_sizes <- c(120)
variable_numbers <- c(100)
R_methods <- c("GHSGEM", "fastGHS")
MATLAB_methods <- c("GHS_MCMC")
structures <- c("random", "bdgraph_sf", "huge_sf", "hubs")

# Read results from the files.
fastGHS_results_p100 <- list()
fastGHS_results_p100 <- add_R_results(fastGHS_results_p100, R_methods, structures, sample_sizes,
                                      variable_numbers)
fastGHS_results_p100 <- add_MATLAB_results(fastGHS_results_p100, MATLAB_methods, structures,
                                           sample_sizes, variable_numbers)

# Set method names.
random_methods <- c("GHSGEM", "GHSGEM_p/2", "GHS_MCMC", "fastGHS")
other_methods <- c("GHSGEM", "GHS_MCMC", "fastGHS")

# Set scores.
scores <- c("MCC", "TPR", "FPR", "FDR", "f_norm_rel", "sl_omega", "time")

# Print results.
print_results(scores, fastGHS_results_p100, random_methods, "random", 120, 100)
print_results(scores, fastGHS_results_p100, other_methods, "bdgraph_sf", 120, 100)
print_results(scores, fastGHS_results_p100, other_methods, "huge_sf", 120, 100)
print_results(scores, fastGHS_results_p100, other_methods, "hubs", 120, 100)

# Results of the fastGHS method with p = 200 and network structure is random.
# Read results from the files.
fastGHS_results_p200 <- list()
fastGHS_results_p200 <- add_R_results(fastGHS_results_p200, R_methods, "random", 120, 200)
fastGHS_results_p200 <- add_MATLAB_results(fastGHS_results_p200, MATLAB_methods, "random", 120, 200)

# Print results.
print_results(scores, fastGHS_results_p200, random_methods, "random", 120, 200)

# Create LaTeX tables.
# Supplementary Table S7.
create_latex_table(scores, fastGHS_results_p100, random_methods, "random", 120, 100)
create_latex_table(scores, fastGHS_results_p100, other_methods, "bdgraph_sf", 120, 100)
create_latex_table(scores, fastGHS_results_p100, other_methods, "huge_sf", 120, 100)
create_latex_table(scores, fastGHS_results_p100, other_methods, "hubs", 120, 100)

# Supplementary Table S8.
create_latex_table(scores, fastGHS_results_p200, random_methods, "random", 120, 200)


######
# Create runtime table. Not presented in the paper or Supplementary Material.
create_runtime_table(all_results, other_methods, structures, sample_sizes, variable_numbers)