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
if(!require("flare", quietly = TRUE)) {
  install.packages("flare")
}
if(!require("GGMncv", quietly = TRUE)) {
  devtools::install_github("donaldRwilliams/GGMncv")
}
if(!require("nlshrink", quietly = TRUE)) {
  install.packages("nlshrink")
}
if(!require("scmamp", quietly = TRUE)) {
  devtools::install_github("b0rxa/scmamp")
}

# Load packages.
library(GHSCM)
library(beam)
library(huge)
library(fastGHS)
library(doParallel)
library(pulsar)
library(flare)
library(GGMncv)
library(nlshrink)

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
R_methods <- c("GHSCM", "GLASSO", "fastGHS", "TIGER", "CLIME", "GSCAD")
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

# Number of empty network estimates of CLIME.
sum(all_results$CLIME$random$n120_p100$results$edge_count == 0)
sum(all_results$CLIME$random$n120_p200$results$edge_count == 0)
sum(all_results$CLIME$bdgraph_sf$n120_p100$results$edge_count == 0)
sum(all_results$CLIME$bdgraph_sf$n120_p200$results$edge_count == 0)

# Set NaNs to zero for CLIME in MCC and FDR columns.
all_results$CLIME$random$n120_p100$results$MCC[is.na(all_results$CLIME$random$n120_p100$results$MCC)] <- 0
all_results$CLIME$random$n120_p100$results$FDR[is.na(all_results$CLIME$random$n120_p100$results$FDR)] <- 0
all_results$CLIME$random$n120_p200$results$MCC[is.na(all_results$CLIME$random$n120_p200$results$MCC)] <- 0
all_results$CLIME$random$n120_p200$results$FDR[is.na(all_results$CLIME$random$n120_p200$results$FDR)] <- 0

all_results$CLIME$bdgraph_sf$n120_p100$results$MCC[is.na(all_results$CLIME$bdgraph_sf$n120_p100$results$MCC)] <- 0
all_results$CLIME$bdgraph_sf$n120_p100$results$FDR[is.na(all_results$CLIME$bdgraph_sf$n120_p100$results$FDR)] <- 0
all_results$CLIME$bdgraph_sf$n120_p200$results$MCC[is.na(all_results$CLIME$bdgraph_sf$n120_p200$results$MCC)] <- 0
all_results$CLIME$bdgraph_sf$n120_p200$results$FDR[is.na(all_results$CLIME$bdgraph_sf$n120_p200$results$FDR)] <- 0

# Define the methods and which order they will be printed. GHS CM has two results for random network
# structure.
random_methods <- c("GHSCM_p/2", "GHSCM", "GHS_MCMC", "fastGHS", "GHS_LLA", "HSL_MCMC", "HSL_ECM",
                    "GLASSO", "TIGER", "CLIME", "GSCAD")
other_methods <- c("GHSCM", "GHS_MCMC", "fastGHS", "GHS_LLA", "HSL_MCMC", "HSL_ECM", "GLASSO",
                   "TIGER", "CLIME", "GSCAD")

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

# Create critical differences plots.
# First combine MCC results into matrices.
column_names <- c("GHS CM", "GHS MCMC", "GHS ECM", "GHS LLA", "GHS-like MCMC", "GHS-like ECM",
                  "GLASSO (StARS)", "TIGER", "CLIME", "GSCAD (BIC)")
mcc_p100 <- mcc_p200 <- matrix(nrow = 50*4, ncol = length(other_methods))
j <- 0
for (structure in structures) {
  for (i in 1:length(other_methods)) {
    mcc_p100[(j*50+1):((j+1)*50),i] <- all_results[[other_methods[i]]][[structure]][["n120_p100"]][["results"]][,"MCC"]
    mcc_p200[(j*50+1):((j+1)*50),i] <- all_results[[other_methods[i]]][[structure]][["n120_p200"]][["results"]][,"MCC"]
  }
  j <- j + 1
}
colnames(results_p100) <- colnames(mcc_p200) <- column_names

# Next define some variables needed for plotting.
text_x_pos <- -0.7
line.spacing <- 0.25
h.up <- 2.5 * line.spacing

# Save figure as eps image.
setEPS()
postscript("Figures/Main_article/CD_diagram.eps", width = 20, height = 5)

par(mfrow = c(1,2))
plotCD(mcc_p100, alpha = 0.05, cex = 1.85, char.size = 0.5)
text(text_x_pos, h.up-0.1, "A", cex = 2.5)
plotCD(mcc_p200, alpha = 0.05, cex = 1.85, char.size = 0.5)
text(text_x_pos, h.up-0.1, "B", cex = 2.5)

dev.off()


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

### Additional visualization of simulation results

methods <- c("GHS CM", "GHS MCMC", "GHS ECM", "GHS LLA", "GHS-like MCMC", "GHS-like ECM",
             "GLASSO (StARS)", "TIGER", "CLIME", "GSCAD (BIC)")

create_boxplots <- function(all_results, p, score, methods, ylim = c(0,1), filename = NULL,
                            image_size = c(10, 7), label_position = -0.03) {
  if (!is.null(filename)) pdf(filename, width = image_size[1], height = image_size[2])
  scores_per_struct <- list()
  for (structure in structures) {
    scores_per_struct[[structure]] <- matrix(nrow = 50, ncol = length(other_methods))
    for (i in 1:length(other_methods)) {
      scores_per_struct[[structure]][,i] <- all_results[[other_methods[i]]][[structure]][[paste0("n120_p", p)]][["results"]][, score]
    }
    colnames(scores_per_struct[[structure]]) <- methods
  }
  par(mfrow = c(2,2))
  par(mar = c(3.1, 8.1, 1.2, 0.2))
  par(mgp = c(2,1,0))
  for (structure in structures) {
    if (structure == "random") title <- "Random network structure"
    else if (structure == "bdgraph_sf") title <- "Scale-free network structure (BDgraph)"
    else if (structure == "huge_sf") title <- "Scale-free network structure (huge)"
    else if (structure == "hubs") title <- "Hub network structure"
    if (score == "f_norm_rel") xlabel <- "Relative F norm"
    else xlabel <- score
    # if (score == "sl_omega" | score ==  "f_norm_omega" | score == "sl_sigma" | score == "f_norm_sigma") {
    #   ylim <- c(min(scores_per_struct[[structure]]), max(scores_per_struct[[structure]]))
    # }
    boxplot(scores_per_struct[[structure]], xaxt = "n", yaxt = "n", xlab = xlabel,
            ylim = ylim, main = title, horizontal = TRUE, at = rev(1:length(methods)))
    
    axis(side = 2, las = 2, labels = FALSE, at = 1:length(methods))
    axis(side = 1)
    
    text(y = 1:length(methods),
         x = par("usr")[1] + label_position * (par("usr")[2] - par("usr")[1]),
         labels = rev(methods),
         xpd = NA,
         cex = 1, adj = 1)
    grid(ny = NA)
  }
  if (!is.null(filename)) dev.off()
}

# Plot only.
create_boxplots(all_results, 100, "MCC", methods, ylim = c(0, 1))
create_boxplots(all_results, 200, "MCC", methods, ylim = c(0, 1))

create_boxplots(all_results, 100, "f_norm_rel", methods, ylim = c(0, 1))
create_boxplots(all_results, 200, "f_norm_rel", methods, ylim = c(0, 1.25))

# Save figures.
create_boxplots(all_results, 100, "MCC", methods, ylim = c(0, 1),
                filename = "Figures/Supplementary/MCC_boxplots_p100.pdf",
                image_size = c(11, 6.5))
create_boxplots(all_results, 200, "MCC", methods, ylim = c(0, 1),
                filename = "Figures/Supplementary/MCC_boxplots_p200.pdf",
                image_size = c(11, 6.5))

create_boxplots(all_results, 100, "f_norm_rel", methods, ylim = c(0, 1),
                filename = "Figures/Supplementary/Rel_f_norm_boxplots_p100.pdf",
                image_size = c(11, 6.5))
create_boxplots(all_results, 200, "f_norm_rel", methods, ylim = c(0, 1.25),
                filename = "Figures/Supplementary/Rel_f_norm_boxplots_p200.pdf",
                label_position = -0.05,
                image_size = c(11, 6.5))

