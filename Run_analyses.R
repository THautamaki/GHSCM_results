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
  else {
    stop("Wrong method!")
  }
  start <- Sys.time()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  results <- foreach (i = 1:n_datasets, .combine = "rbind", .packages = packages) %dopar% {
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
  stopCluster(cl)
  stopImplicitCluster()
  total_time <- Sys.time() - start
  return(list(results = results, total_time = total_time))
}

ps <- c(100, 200)

# Analysis for the datasets with 100 variables.

p <- 100
path <- paste0("C:\\Users\\tume8\\Documents\\data\\Artikkeli\\p", p, "\\")

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

scores <- c("MCC", "F1", "TPR", "FPR", "FDR", "sl_omega", "F_norm_omega", "time")

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

save(random_p100_results, file = "random_p100_GEM_results.Rda")
save(bdgraph_sf_p100_results, file = "bdgraph_sf_p100_GEM_results.Rda")
save(huge_sf_p100_results, file = "huge_sf_p100_GEM_results.Rda")
save(hubs_p100_results, file = "hubs_p100_GEM_results.Rda")

save(random_p100_glasso_results, file = "random_p100_glasso_results.Rda")
save(bdgraph_sf_p100_glasso_results, file = "bdgraph_sf_p100_glasso_results.Rda")
save(huge_sf_p100_glasso_results, file = "huge_sf_p100_glasso_results.Rda")
save(hubs_p100_glasso_results, file = "hubs_p100_glasso_results.Rda")

save(random_p100_beam_results, file = "random_p100_beam_results.Rda")
save(bdgraph_sf_p100_beam_results, file = "bdgraph_sf_p100_beam_results.Rda")
save(huge_sf_p100_beam_results, file = "huge_sf_p100_beam_results.Rda")
save(hubs_p100_beam_results, file = "hubs_p100_beam_results.Rda")

# Analysis for the datasets with 200 variables.

p <- 200
path <- paste0("C:\\Users\\tume8\\Documents\\data\\Artikkeli\\p", p, "\\")

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

save(random_p200_results, file = "random_p200_GEM_results.Rda")
save(bdgraph_sf_p200_results, file = "bdgraph_sf_p200_GEM_results.Rda")
save(huge_sf_p200_results, file = "huge_sf_p200_GEM_results.Rda")
save(hubs_p200_results, file = "hubs_p200_GEM_results.Rda")

save(random_p200_glasso_results, file = "random_p200_glasso_results.Rda")
save(bdgraph_sf_p200_glasso_results, file = "bdgraph_sf_p200_glasso_results.Rda")
save(huge_sf_p200_glasso_results, file = "huge_sf_p200_glasso_results.Rda")
save(hubs_p200_glasso_results, file = "hubs_p200_glasso_results.Rda")

save(random_p200_beam_results, file = "random_p200_beam_results.Rda")
save(bdgraph_sf_p200_beam_results, file = "bdgraph_sf_p200_beam_results.Rda")
save(huge_sf_p200_beam_results, file = "huge_sf_p200_beam_results.Rda")
save(hubs_p200_beam_results, file = "hubs_p200_beam_results.Rda")


########

n_datasets <- length(sim_random)
cores <- detectCores(logical = FALSE)

start <- Sys.time()
cl <- makeCluster(cores)
registerDoParallel(cl)
results <- foreach (i = 1:n_datasets, .combine = "rbind", .packages = "GHSGEM") %dopar% {
  sim <- sim_random[[i]]
  map <- GHS_MAP_estimate(sim$data, verbose = 0)
  theta_est <- map$Theta
  sigma_est <- map$Sigma_est
  omega_est <- map$Omega_est
  time <- map$tot_time

  edge_count <- sum(theta_est) / 2
  cm <- conf_matrix(sim$theta, theta_est)
  sl_omega <- stein_loss(sim$omega, omega_est)
  sl_sigma <- stein_loss(sim$sigma, sigma_est)
  F_norm_sigma <- norm(sim$sigma - sigma_est, type = "f")
  F_norm_omega <- norm(sim$omega - omega_est, type = "f")
  score <- cbind(i, calculate_scores(cm), edge_count, sl_sigma, sl_omega, F_norm_sigma,
                 F_norm_omega, time)
  score
}
stopCluster(cl)
stopImplicitCluster()
Sys.time() - start

