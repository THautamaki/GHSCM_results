if(!require("GHSGEM", quietly = TRUE)) {
  devtools::install_github("THautamaki/GHSGEM")
}
if(!require("beam", quietly = TRUE)) {
  devtools::install_version("beam", version = "2.0.2")
}

library(GHSGEM)
library(beam)
library(huge)
library(doParallel)
library(pulsar)

run_analysis <- function(simulations, method = "GEM") {
  n_datasets <- length(simulations)
  cores <- detectCores(logical = FALSE)
  if (tolower(method) == "gem") {
    packages <- c("GHSGEM")
  }
  else if (tolower(method) == "glasso") {
    packages <- c("huge", "GHSGEM")
  }
  else if (tolower(method) == "beam") {
    packages <- c("beam", "GHSGEM")
  }
  else if (tolower(method) == "pulsar") {
    packages <- c("pulsar", "huge", "GHSGEM")
  }
  else {
    stop("Wrong method!")
  }
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  start <- Sys.time()
  results <- foreach (i = 1:n_datasets, .combine = "rbind", .packages = packages,
                      .verbose = FALSE) %dopar% {
    sim <- simulations[[i]]
    if (tolower(method) == "gem") {
      map <- GHS_MAP_estimate(sim$data, verbose = 0)
      theta_est <- map$Theta
      sigma_est <- map$Sigma_est
      omega_est <- map$Omega_est
      time <- map$tot_time
    }
    else if (tolower(method) == "glasso") {
      start_glasso <- Sys.time()
      glasso_est <- huge(sim$data, nlambda = 50, method = "glasso")
      stars_select <- huge.select(glasso_est, criterion = "stars")
      time <- as.numeric(difftime(Sys.time(), start_glasso, units = "secs"))
      theta_est <- stars_select$refit
      omega_est <- stars_select$opt.icov
    }
    else if (tolower(method) == "beam") {
      beam_start <- Sys.time()
      beam_est <- beam(sim$data)
      beam_select <- beam.select(beam_est)
      time <- as.numeric(difftime(Sys.time(), beam_start, units = "secs"))
      theta_est <- igraph::as_adjacency_matrix(ugraph(beam_select), names = FALSE)
    }
    else if (tolower(method) == "pulsar") {
      pulsar_start <- Sys.time()
      #lmax <- getMaxCov(cov(sim$data))
      #lams <- getLamPath(lmax, lmax * 0.05, len = 50)
      lams <- huge(sim$data, nlambda = 50, method = "glasso")$lambda
      bdgraphargs <- list(lambda = lams, method = "glasso", verbose = FALSE)
      out.p <- pulsar(sim$data, fun = huge, fargs = bdgraphargs, rep.num = 20, criterion = 'stars',
                      lb.stars = TRUE, ub.stars = TRUE)
      fit.p <- refit(out.p)
      time <- as.numeric(difftime(Sys.time(), pulsar_start, units = "secs"))
      theta_est <- as.matrix(fit.p$refit$stars)
      diag(theta_est) <- 0
    }
    edge_count <- sum(theta_est) / 2
    cm <- conf_matrix(sim$theta, theta_est)
    sl_omega <- NA
    sl_sigma <- NA
    F_norm_sigma <- NA
    F_norm_omega <- NA
    if (tolower(method) == "gem") {
      sl_omega <- stein_loss(sim$omega, omega_est)
      sl_sigma <- stein_loss(sim$sigma, sigma_est)
      F_norm_omega <- norm(sim$omega - omega_est, type = "f")
      F_norm_sigma <- norm(sim$sigma - sigma_est, type = "f")
    }
    else if (tolower(method) == "glasso") {
      sl_omega <- stein_loss(sim$omega, omega_est)
      F_norm_omega <- norm(sim$omega - omega_est, type = "f")
    }
    score <- cbind(i, calculate_scores(cm), edge_count, sl_sigma, sl_omega, F_norm_sigma,
                   F_norm_omega, time)
    score
  }
  total_time <- Sys.time() - start
  stopCluster(cl)
  stopImplicitCluster()
  return(list(results = results, total_time = total_time))
}

# Analysis for the datasets with 100 variables.
p <- 100

# Change "path\\to\\data" part where you have stored datasets.
path <- paste0("path\\to\\data\\p", p, "\\")

load(file = paste0(path, "bdgraph_random_", p, ".Rda"))
load(file = paste0(path, "bdgraph_scale-free_", p, ".Rda"))
load(file = paste0(path, "huge_scale-free_", p, ".Rda"))
load(file = paste0(path, "huge_hubs_", p, ".Rda"))

random_p100_results <- run_analysis(sim_random)
bdgraph_sf_p100_results <- run_analysis(sim_bdgraph_sf)
huge_sf_p100_results <- run_analysis(sim_huge_sf)
hubs_p100_results <- run_analysis(sim_hubs)

random_p100_beam_results <- run_analysis(sim_random, "beam")
bdgraph_sf_p100_beam_results <- run_analysis(sim_bdgraph_sf, "beam")
huge_sf_p100_beam_results <- run_analysis(sim_huge_sf, "beam")
hubs_p100_beam_results <- run_analysis(sim_hubs, "beam")

random_p100_glasso_results <- run_analysis(sim_random, "glasso")
bdgraph_sf_p100_glasso_results <- run_analysis(sim_bdgraph_sf, "glasso")
huge_sf_p100_glasso_results <- run_analysis(sim_huge_sf, "glasso")
hubs_p100_glasso_results <- run_analysis(sim_hubs, "glasso")

bdgraph_sf_p100_pulsar_results <- run_analysis(sim_bdgraph_sf, "pulsar")
huge_sf_p100_pulsar_results <- run_analysis(sim_huge_sf, "pulsar")

# Select scores which wanted to print.
scores <- c("MCC", "F1", "TPR", "FPR", "FDR", "sl_omega", "F_norm_omega", "time")

# Print scores.
round(rbind(colMeans(random_p100_results$results[, scores]),
            colMeans(random_p100_glasso_results$results[, scores]),
            colMeans(random_p100_beam_results$results[, scores])), 4)

round(rbind(colMeans(bdgraph_sf_p100_results$results[, scores]),
            colMeans(bdgraph_sf_p100_glasso_results$results[, scores]),
            colMeans(bdgraph_sf_p100_beam_results$results[, scores])), 4)

round(rbind(colMeans(huge_sf_p100_results$results[, scores]),
            colMeans(huge_sf_p100_glasso_results$results[, scores]),
            colMeans(huge_sf_p100_beam_results$results[, scores])), 4)

round(rbind(colMeans(hubs_p100_results$results[, scores]),
            colMeans(hubs_p100_glasso_results$results[, scores]),
            colMeans(hubs_p100_beam_results$results[, scores])), 4)

# Save results into "Results_files" folder which is in the working directory.
save(random_p100_results, file = "Results_files/random_p100_GEM_results.Rda")
save(bdgraph_sf_p100_results, file = "Results_files/bdgraph_sf_p100_GEM_results.Rda")
save(huge_sf_p100_results, file = "Results_files/huge_sf_p100_GEM_results.Rda")
save(hubs_p100_results, file = "Results_files/hubs_p100_GEM_results.Rda")

save(random_p100_glasso_results, file = "Results_files/random_p100_glasso_results.Rda")
save(bdgraph_sf_p100_glasso_results, file = "Results_files/bdgraph_sf_p100_glasso_results.Rda")
save(huge_sf_p100_glasso_results, file = "Results_files/huge_sf_p100_glasso_results.Rda")
save(hubs_p100_glasso_results, file = "Results_files/hubs_p100_glasso_results.Rda")

save(random_p100_beam_results, file = "Results_files/random_p100_beam_results.Rda")
save(bdgraph_sf_p100_beam_results, file = "Results_files/bdgraph_sf_p100_beam_results.Rda")
save(huge_sf_p100_beam_results, file = "Results_files/huge_sf_p100_beam_results.Rda")
save(hubs_p100_beam_results, file = "Results_files/hubs_p100_beam_results.Rda")


# Analysis for the datasets with 200 variables.
p <- 200

# Change "path\\to\\data" part where you have stored datasets.
path <- paste0("path\\to\\data\\p", p, "\\")

load(file = paste0(path, "bdgraph_random_", p, ".Rda"))
load(file = paste0(path, "bdgraph_scale-free_", p, ".Rda"))
load(file = paste0(path, "huge_scale-free_", p, ".Rda"))
load(file = paste0(path, "huge_hubs_", p, ".Rda"))

n_datasets <- length(sim_random)

random_p200_results <- run_analysis(sim_random)
bdgraph_sf_p200_results <- run_analysis(sim_bdgraph_sf)
huge_sf_p200_results <- run_analysis(sim_huge_sf)
hubs_p200_results <- run_analysis(sim_hubs)

random_p200_beam_results <- run_analysis(sim_random, "beam")
bdgraph_sf_p200_beam_results <- run_analysis(sim_bdgraph_sf, "beam")
huge_sf_p200_beam_results <- run_analysis(sim_huge_sf, "beam")
hubs_p200_beam_results <- run_analysis(sim_hubs, "beam")

random_p200_glasso_results <- run_analysis(sim_random, "glasso")
bdgraph_sf_p200_glasso_results <- run_analysis(sim_bdgraph_sf, "glasso")
huge_sf_p200_glasso_results <- run_analysis(sim_huge_sf, "glasso")
hubs_p200_glasso_results <- run_analysis(sim_hubs, "glasso")

# Select scores which wanted to print.
scores <- c("MCC", "F1", "TPR", "FPR", "FDR", "sl_omega", "F_norm_omega", "time")

# Print scores.
round(rbind(colMeans(random_p200_results$results[, scores]),
            colMeans(random_p200_glasso_results$results[, scores]),
            colMeans(random_p200_beam_results$results[, scores])), 4)

round(rbind(colMeans(bdgraph_sf_p200_results$results[, scores]),
            colMeans(bdgraph_sf_p200_glasso_results$results[, scores]),
            colMeans(bdgraph_sf_p200_beam_results$results[, scores])), 4)

round(rbind(colMeans(huge_sf_p200_results$results[, scores]),
            colMeans(huge_sf_p200_glasso_results$results[, scores]),
            colMeans(huge_sf_p200_beam_results$results[, scores])), 4)

round(rbind(colMeans(hubs_p200_results$results[, scores]),
            colMeans(hubs_p200_glasso_results$results[, scores]),
            colMeans(hubs_p200_beam_results$results[, scores])), 4)

save(random_p200_results, file = "Results_files/random_p200_GEM_results.Rda")
save(bdgraph_sf_p200_results, file = "Results_files/bdgraph_sf_p200_GEM_results.Rda")
save(huge_sf_p200_results, file = "Results_files/huge_sf_p200_GEM_results.Rda")
save(hubs_p200_results, file = "Results_files/hubs_p200_GEM_results.Rda")

save(random_p200_glasso_results, file = "Results_files/random_p200_glasso_results.Rda")
save(bdgraph_sf_p200_glasso_results, file = "Results_files/bdgraph_sf_p200_glasso_results.Rda")
save(huge_sf_p200_glasso_results, file = "Results_files/huge_sf_p200_glasso_results.Rda")
save(hubs_p200_glasso_results, file = "Results_files/hubs_p200_glasso_results.Rda")

save(random_p200_beam_results, file = "Results_files/random_p200_beam_results.Rda")
save(bdgraph_sf_p200_beam_results, file = "Results_files/bdgraph_sf_p200_beam_results.Rda")
save(huge_sf_p200_beam_results, file = "Results_files/huge_sf_p200_beam_results.Rda")
save(hubs_p200_beam_results, file = "Results_files/hubs_p200_beam_results.Rda")


#########
# Recalculate scores if analyses are already run and saved.
# Load results files.

load(file = "Results_files/random_p100_GEM_results.Rda")
load(file = "Results_files/bdgraph_sf_p100_GEM_results.Rda")
load(file = "Results_files/huge_sf_p100_GEM_results.Rda")
load(file = "Results_files/hubs_p100_GEM_results.Rda")

load(file = "Results_files/random_p100_glasso_results.Rda")
load(file = "Results_files/bdgraph_sf_p100_glasso_results.Rda")
load(file = "Results_files/huge_sf_p100_glasso_results.Rda")
load(file = "Results_files/hubs_p100_glasso_results.Rda")

load(file = "Results_files/random_p100_beam_results.Rda")
load(file = "Results_files/bdgraph_sf_p100_beam_results.Rda")
load(file = "Results_files/huge_sf_p100_beam_results.Rda")
load(file = "Results_files/hubs_p100_beam_results.Rda")

load(file = "Results_files/random_p200_GEM_results.Rda")
load(file = "Results_files/bdgraph_sf_p200_GEM_results.Rda")
load(file = "Results_files/huge_sf_p200_GEM_results.Rda")
load(file = "Results_files/hubs_p200_GEM_results.Rda")

load(file = "Results_files/random_p200_glasso_results.Rda")
load(file = "Results_files/bdgraph_sf_p200_glasso_results.Rda")
load(file = "Results_files/huge_sf_p200_glasso_results.Rda")
load(file = "Results_files/hubs_p200_glasso_results.Rda")

load(file = "Results_files/random_p200_beam_results.Rda")
load(file = "Results_files/bdgraph_sf_p200_beam_results.Rda")
load(file = "Results_files/huge_sf_p200_beam_results.Rda")
load(file = "Results_files/hubs_p200_beam_results.Rda")

scores <- c("MCC", "F1", "TPR", "FPR", "FDR", "sl_omega", "F_norm_omega", "time")

# Means of the scores:
# p = 100.
# Random network.
round(rbind(colMeans(random_p100_results$results[, scores]),
            colMeans(random_p100_glasso_results$results[, scores]),
            colMeans(random_p100_beam_results$results[, scores])), 4)

# Scale-free network (BDgraph).
round(rbind(colMeans(bdgraph_sf_p100_results$results[, scores]),
            colMeans(bdgraph_sf_p100_glasso_results$results[, scores]),
            colMeans(bdgraph_sf_p100_beam_results$results[, scores])), 4)

# Scale-free network (huge).
round(rbind(colMeans(huge_sf_p100_results$results[, scores]),
            colMeans(huge_sf_p100_glasso_results$results[, scores]),
            colMeans(huge_sf_p100_beam_results$results[, scores])), 4)

# Hubs network.
round(rbind(colMeans(hubs_p100_results$results[, scores]),
            colMeans(hubs_p100_glasso_results$results[, scores]),
            colMeans(hubs_p100_beam_results$results[, scores])), 4)

# p = 200.
# Random network.
round(rbind(colMeans(random_p200_results$results[, scores]),
            colMeans(random_p200_glasso_results$results[, scores]),
            colMeans(random_p200_beam_results$results[, scores])), 4)

# Scale-free network (BDgraph).
round(rbind(colMeans(bdgraph_sf_p200_results$results[, scores]),
            colMeans(bdgraph_sf_p200_glasso_results$results[, scores]),
            colMeans(bdgraph_sf_p200_beam_results$results[, scores])), 4)

# Scale-free network (huge).
round(rbind(colMeans(huge_sf_p200_results$results[, scores]),
            colMeans(huge_sf_p200_glasso_results$results[, scores]),
            colMeans(huge_sf_p200_beam_results$results[, scores])), 4)

# Hubs network.
round(rbind(colMeans(hubs_p200_results$results[, scores]),
            colMeans(hubs_p200_glasso_results$results[, scores]),
            colMeans(hubs_p200_beam_results$results[, scores])), 4)

# SD:s of the results:
# p = 100.
# Random network.
round(rbind(apply(random_p100_results$results[, scores], 2, sd),
            apply(random_p100_glasso_results$results[, scores], 2, sd),
            apply(random_p100_beam_results$results[, scores], 2, sd)), 4)

# Scale-free network (BDgraph).
round(rbind(apply(bdgraph_sf_p100_results$results[, scores], 2, sd),
            apply(bdgraph_sf_p100_glasso_results$results[, scores], 2, sd),
            apply(bdgraph_sf_p100_beam_results$results[, scores], 2, sd)), 4)

# Scale-free network (huge).
round(rbind(apply(huge_sf_p100_results$results[, scores], 2, sd),
            apply(huge_sf_p100_glasso_results$results[, scores], 2, sd),
            apply(huge_sf_p100_beam_results$results[, scores], 2, sd)), 4)

# Hubs network.
round(rbind(apply(hubs_p100_results$results[, scores], 2, sd),
            apply(hubs_p100_glasso_results$results[, scores], 2, sd),
            apply(hubs_p100_beam_results$results[, scores], 2, sd)), 4)

# p = 200.
# Random network.
round(rbind(apply(random_p200_results$results[, scores], 2, sd),
            apply(random_p200_glasso_results$results[, scores], 2, sd),
            apply(random_p200_beam_results$results[, scores], 2, sd)), 4)

# Scale-free network (BDgraph).
round(rbind(apply(bdgraph_sf_p200_results$results[, scores], 2, sd),
            apply(bdgraph_sf_p200_glasso_results$results[, scores], 2, sd),
            apply(bdgraph_sf_p200_beam_results$results[, scores], 2, sd)), 4)

# Scale-free network (huge).
round(rbind(apply(huge_sf_p200_results$results[, scores], 2, sd),
            apply(huge_sf_p200_glasso_results$results[, scores], 2, sd),
            apply(huge_sf_p200_beam_results$results[, scores], 2, sd)), 4)

# Hubs network.
round(rbind(apply(hubs_p200_results$results[, scores], 2, sd),
            apply(hubs_p200_glasso_results$results[, scores], 2, sd),
            apply(hubs_p200_beam_results$results[, scores], 2, sd)), 4)

# Total running times:
# p = 100.
# Random network.
random_p100_results$total_time
random_p100_glasso_results$total_time
random_p100_beam_results$total_time

# Scale-free network (BDgraph).
bdgraph_sf_p100_results$total_time
as.double(bdgraph_sf_p100_glasso_results$total_time, units = "secs")
bdgraph_sf_p100_beam_results$total_time

# Scale-free network (huge).
huge_sf_p100_results$total_time
huge_sf_p100_glasso_results$total_time
huge_sf_p100_beam_results$total_time

# Hubs network.
hubs_p100_results$total_time
hubs_p100_glasso_results$total_time
hubs_p100_beam_results$total_time

# p = 200.
# Random network.
random_p200_results$total_time
as.double(random_p200_glasso_results$total_time, units = "secs")
random_p200_beam_results$total_time

# Scale-free network (BDgraph).
as.double(bdgraph_sf_p200_results$total_time, units = "secs")
as.double(bdgraph_sf_p200_glasso_results$total_time, units = "secs")
bdgraph_sf_p200_beam_results$total_time

# Scale-free network (huge).
as.double(huge_sf_p200_results$total_time, units = "secs")
as.double(huge_sf_p200_glasso_results$total_time, units = "secs")
huge_sf_p200_beam_results$total_time

# Hubs network.
as.double(hubs_p200_results$total_time, units = "secs")
as.double(hubs_p200_glasso_results$total_time, units = "secs")
hubs_p200_beam_results$total_time

