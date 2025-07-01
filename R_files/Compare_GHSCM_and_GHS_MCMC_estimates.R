library(GHSGEM)
library(R.matlab)
library(doParallel)

compare_GHSCM_and_GHS_MCMC <- function(structures, sample_sizes, variable_numbers, save_results = TRUE) {
  run_analysis <- function(structure, n, p, p0 = 0, save_results = TRUE) {
    path <- paste0("Data/n", n, "_p", p, "/")
    if (tolower(structure) == "random") {
      simulations <- readRDS(file = paste0(path, "bdgraph_random_n", n, "_p", p, ".Rds"))
    }
    else if (tolower(structure) == "bdgraph_sf") {
      simulations <- readRDS(file = paste0(path, "bdgraph_scale-free_n", n, "_p", p, ".Rds"))
    }
    else if (tolower(structure) == "huge_sf") {
      simulations <- readRDS(file = paste0(path, "huge_scale-free_n", n, "_p", p, ".Rds"))
    }
    else if (tolower(structure) == "hubs") {
      simulations <- readRDS(file = paste0(path, "huge_hubs_n", n, "_p", p, ".Rds"))
    }
    else {
      stop("Wrong network structure! Possible choices are 'random', 'bdgraph_sf', 'huge_sf' and 'hubs'.")
    }
    n_datasets <- length(simulations)
    cl <- makeCluster(detectCores(logical = FALSE))
    registerDoParallel(cl)
    maps <- foreach (r = 1:n_datasets, .packages = "GHSGEM") %dopar% {
      GHS_MAP_estimation(simulations[[r]]$data, p0 = p0)
    }
    stopCluster(cl)
    stopImplicitCluster()
    path_matrices <- paste0("Results_files/p", p, "/GHS_MCMC/", structure,  "/")
    results <- data.frame()
    for (r in 1:n_datasets) {
      matrices <- readMat(paste0(path_matrices, structure,  "_matrices_", r, ".mat"))
      cm <- conf_matrix(matrices$Theta, maps[[r]]$Theta_est)
      sl_mean <- stein_loss(matrices$Omega.mean, maps[[r]]$Omega_est)
      f_norm_mean <- norm(matrices$Omega.mean - maps[[r]]$Omega_est, type = "F")
      sl_median <- stein_loss(matrices$Omega.median, maps[[r]]$Omega_est)
      f_norm_median <- norm(matrices$Omega.median - maps[[r]]$Omega_est, type = "F")
      results <- rbind(results, cbind(r, calculate_scores(cm), sl_mean, f_norm_mean, sl_median, f_norm_median))
    }
    results <- list(results = results)
    if (save_results) {
      results_path <- paste0("Results_files/p", p, "/GHSCM_vs_GHS_MCMC/")
      if (tolower(structure) == "random" & p0 > 0) {
        saveRDS(results, file = paste0(results_path, structure, "_p", p, "_GHSCM_vs_GHS_MCMC_results_2.Rds"))
      }
      else {
        saveRDS(results, file = paste0(results_path, structure, "_p", p, "_GHSCM_vs_GHS_MCMC_results.Rds"))
      }
    }
    return(results)
  }
  all_results <- list()
  for (structure in structures) {
    for (n in sample_sizes) {
      for (p in variable_numbers) {
        all_results[[structure]][["GHSCM_vs_GHS_MCMC"]][[paste0("n", n, "_p", p)]] <- run_analysis(structure, n, p, save_results = save_results)
        if (structure == "random") {
          all_results[[paste0(structure, "_p/2")]][["GHSCM_vs_GHS_MCMC"]][[paste0("n", n, "_p", p)]] <- run_analysis(structure, n, p, p0 = p/2, save_results = save_results)
        }
      }
    }
  }
  return(all_results)
}

add_R_results2 <- function(results_list, methods, structures, sample_sizes, variable_numbers) {
  for (method in methods) {
    for (structure in structures) {
      for (n in sample_sizes) {
        for (p in variable_numbers) {
          path <- paste0("Results_files/p", p, "/", structure, "/")
          if (tolower(method) == "random") {
            results <- readRDS(paste0(path, method, "_p", p, "_", structure, "_results_2.Rds"))
            results_list[[paste0(method, "_p/2")]][[structure]][[paste0("n", n, "_p", p)]]$results <- results$results
          }
          path <- paste0("Results_files/p", p, "/", structure, "/")
          results <- readRDS(paste0(path, method, "_p", p, "_", structure, "_results.Rds"))
          results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results <- results$results
        }
      }
    }
  }
  return(results_list)
}

# Set all needed parameters.
sample_sizes <- c(120)
variable_numbers <- c(100, 200)
structures <- c("random", "bdgraph_sf", "huge_sf", "hubs")

# Run analysis if not were run.
GHSCM_vs_GHS_MCMC_results <- compare_GHSCM_and_GHS_MCMC(structures = structures, sample_sizes = sample_sizes,
                                                        variable_numbers = variable_numbers, save_results = TRUE)

# Read results from the files, if analysis already done.
GHSCM_vs_GHS_MCMC_results <- list()
GHSCM_vs_GHS_MCMC_results <- add_R_results2(GHSCM_vs_GHS_MCMC_results, methods = structures,
                                            structures = "GHSCM_vs_GHS_MCMC", sample_sizes = sample_sizes,
                                            variable_numbers = variable_numbers)

# Define scores.
scores <- c("MCC", "TPR", "FDR", "sl_mean")
# Print latex tables (Table 3 in the paper).
create_latex_table(scores, GHSCM_vs_GHS_MCMC_results, c("random_p/2", structures),
                   structure = "GHSCM_vs_GHS_MCMC", n = 120, p = 100, highlight_best = FALSE)
create_latex_table(scores, GHSCM_vs_GHS_MCMC_results, c("random_p/2", structures),
                   structure = "GHSCM_vs_GHS_MCMC", n = 120, p = 200, highlight_best = FALSE)
