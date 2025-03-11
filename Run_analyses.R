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

run_analysis <- function(structure, n, p, method = "GEM", p0 = 0) {
  # Load datasets from file.
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
  else if (tolower(method) == "fastghs") {
    packages <- c("fastGHS", "GHSGEM")
  }
  else {
    stop("Wrong method! Possible choices are 'GEM', 'GLASSO', 'pulsar' and 'fastGHS'.")
  }
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  start <- Sys.time()
  results <- foreach (i = 1:n_datasets, .combine = "rbind", .packages = packages,
                      .verbose = FALSE) %dopar% {
    sim <- simulations[[i]]
    if (tolower(method) == "gem") {
      map <- GHS_MAP_estimation(sim$data, verbose = 0, p0 = p0, max_iterations = 1000, tol = 1e-4)
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
    else if (tolower(method) == "fastghs") {
      n <- nrow(sim$data)
      p <- ncol(sim$data)
      fastghs_start <- Sys.time()
      result <- fastGHS(sim$data, AIC_selection = TRUE, epsilon = 1e-3)
      time <- as.numeric(difftime(Sys.time(), fastghs_start, units = "secs"))
      theta_est <- result$theta
      theta_est[abs(theta_est) < 1e-5] <- 0
      theta_est[theta_est != 0] <- 1
      diag(theta_est) <- 0
      omega_est <- result$theta
      sigma_est <- result$sigma
    }
    edge_count <- sum(theta_est) / 2
    cm <- conf_matrix(sim$theta, theta_est)
    sl_omega <- NA
    sl_sigma <- NA
    F_norm_sigma <- NA
    F_norm_omega <- NA
    F_norm_rel <- NA
    if (tolower(method) == "gem" | tolower(method) == "fastghs") {
      sl_omega <- stein_loss(sim$omega, omega_est)
      sl_sigma <- stein_loss(sim$sigma, sigma_est)
      F_norm_omega <- norm(sim$omega - omega_est, type = "f")
      F_norm_sigma <- norm(sim$sigma - sigma_est, type = "f")
      F_norm_rel <- F_norm_omega / norm(sim$omega, type = "f")
    }
    else if (tolower(method) == "glasso") {
      sl_omega <- stein_loss(sim$omega, omega_est)
      F_norm_omega <- norm(sim$omega - omega_est, type = "f")
    }
    score <- cbind(i, calculate_scores(cm), edge_count, sl_sigma, sl_omega, F_norm_sigma,
                   F_norm_omega, F_norm_rel, time)
    score
  }
  total_time <- Sys.time() - start
  stopCluster(cl)
  stopImplicitCluster()
  return(list(results = results, total_time = total_time))
}

structures <- c("random", "bdgraph_sf", "huge_sf", "hubs")

# Analysis for the datasets with 100 variables.
n <- 120
p <- 200

# Change "path\\to\\data" part where you have stored datasets.
path <- paste0("path\\to\\data\\p", p, "\\")

path <- "C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\data\\Artikkeli\\p100\\"
path <- paste0("Data/n", n, "_p", p, "/")

load(file = paste0(path, "/bdgraph_random_n", n, "_p", p, ".Rda"))
load(file = paste0(path, "/bdgraph_scale-free_", p, ".Rda"))
load(file = paste0(path, "/huge_scale-free_", p, ".Rda"))
load(file = paste0(path, "/huge_hubs_", p, ".Rda"))

random_p100_results <- run_analysis(sim_random)
random_p100_results_2 <- run_analysis(sim_random, p0 = p/2)
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

p0 <- p * (p - 1) / 2
n <- 120
p0 <- p*100
p0 / (p * (p-1) / 2) * (50 * sqrt(p0) / (n * sqrt(n)))

random_p100_fastghs_results <- run_analysis(sim_random, "fastGHS", p0 = p*2)
bdgraph_sf_p100_fastghs_results <- run_analysis(sim_bdgraph_sf, "fastGHS", p0 = p*2)
huge_sf_p100_fastghs_results <- run_analysis(sim_huge_sf, "fastGHS", p0 = p*100)
hubs_p100_fastghs_results <- run_analysis(sim_hubs, "fastGHS", p0 = p*30)

# Select scores which will be printed.
scores <- c("MCC", "TPR", "FPR", "FDR", "sl_omega", "F_norm_omega", "F_norm_omega_normalised", "time")
scores <- c("MCC", "sl_omega")

round(rbind(colMeans(random_p100_results$results[, scores]),
            colMeans(random_p100_results_2$results[, scores]),
            colMeans(GHS_MCMC_p100_random[, scores]),
            colMeans(random_p100_fastghs_results$results[,scores])), 4)

round(rbind(colMeans(bdgraph_sf_p100_results$results[, scores]),
            colMeans(GHS_MCMC_p100_bdgraph_sf[, scores])), 4)

round(rbind(colMeans(huge_sf_p100_results$results[, scores]),
            colMeans(GHS_MCMC_p100_huge_sf[, scores])), 4)

round(rbind(colMeans(hubs_p100_results$results[, scores]),
            colMeans(GHS_MCMC_p100_hubs[, scores])), 4)

round(colMeans(random_p100_results$results[, scores]), 4)
round(apply(random_p100_results$results[, scores], 2, sd), 4)

round(colMeans(random_p100_results$results[, scores]), 4)
round(colMeans(random_p100_results_2$results[, scores]), 4)
round(colMeans(random_p100_fastghs_results$results[, scores]), 4)

round(colMeans(bdgraph_sf_p100_results$results[, scores]), 4)
round(colMeans(bdgraph_sf_p100_fastghs_results$results[, scores]), 4)

round(colMeans(huge_sf_p100_results$results[, scores]), 4)
round(colMeans(huge_sf_p100_fastghs_results$results[, scores]), 4)

round(colMeans(hubs_p100_results$results[, scores]), 4)
round(colMeans(hubs_p100_fastghs_results$results[, scores]), 4)

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


# Create table for article usage.
GEM_p100_means <- rbind(colMeans(random_p100_results$results[, scores]),
                        colMeans(bdgraph_sf_p100_results$results[, scores]),
                        colMeans(huge_sf_p100_results$results[, scores]),
                        colMeans(hubs_p100_results$results[, scores]))

GEM_p100_sds <- rbind(apply(random_p100_results$results[, scores], 2, sd),
                      apply(bdgraph_sf_p100_results$results[, scores], 2, sd),
                      apply(huge_sf_p100_results$results[, scores], 2, sd),
                      apply(hubs_p100_results$results[, scores], 2, sd))

GEM_p100_results <- data.frame(row.names = structures)

for (i in 1:ncol(GEM_p100_means)) {
  GEM_p100_results <- cbind(GEM_p100_results, GEM_p100_means[,i], GEM_p100_sds[,i])
}
columns <- c("MCC", "sd", "TPR", "sd", "FPR", "sd", "FDR", "sd", "sl_omega", "sd", "f_norm", "sd", "time", "sd")
colnames(GEM_p100_results) <- columns
round(GEM_p100_results, 3)

round(GEM_p100_results[, 5:6], 4)


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

random_p200_results <- run_analysis(sim_random, p0 = p/2)
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

######
scores <- c("MCC", "TPR", "FPR", "FDR", "sl_omega", "F_norm_omega", "F_norm_omega_normalised", "time")

round(colMeans(random_p200_results$results[, scores]), 4)
round(apply(random_p200_results$results[, scores], 2, sd), 4)

# Select scores which wanted to print.
scores <- c("MCC", "TPR", "FPR", "FDR", "sl_omega", "F_norm_omega", "time")

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

# Create table for article usage.
GEM_p200_means <- rbind(colMeans(random_p200_results$results[, scores]),
                        colMeans(bdgraph_sf_p200_results$results[, scores]),
                        colMeans(huge_sf_p200_results$results[, scores]),
                        colMeans(hubs_p200_results$results[, scores]))

GEM_p200_sds <- rbind(apply(random_p200_results$results[, scores], 2, sd),
                      apply(bdgraph_sf_p200_results$results[, scores], 2, sd),
                      apply(huge_sf_p200_results$results[, scores], 2, sd),
                      apply(hubs_p200_results$results[, scores], 2, sd))

GEM_p200_results <- data.frame(row.names = structures)

for (i in 1:ncol(GEM_p100_means)) {
  GEM_p200_results <- cbind(GEM_p200_results, GEM_p200_means[,i], GEM_p200_sds[,i])
}
columns <- c("MCC", "sd", "TPR", "sd", "FPR", "sd", "FDR", "sd", "sl_omega", "sd", "f_norm", "sd", "time", "sd")
colnames(GEM_p200_results) <- columns
round(GEM_p200_results, 3)

round(GEM_p200_results[, 5:6], 4)

######

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

load(file = "Results_files/p100/GHSGEM/random_p100_GEM_results.Rda")
load(file = "Results_files/p100/GHSGEM/bdgraph_sf_p100_GEM_results.Rda")
load(file = "Results_files/p100/GHSGEM/huge_sf_p100_GEM_results.Rda")
load(file = "Results_files/p100/GHSGEM/hubs_p100_GEM_results.Rda")

load(file = "Results_files/p100/glasso/random_p100_glasso_results.Rda")
load(file = "Results_files/p100/glasso/bdgraph_sf_p100_glasso_results.Rda")
load(file = "Results_files/p100/glasso/huge_sf_p100_glasso_results.Rda")
load(file = "Results_files/p100/glasso/hubs_p100_glasso_results.Rda")

load(file = "Results_files/p100/beam/random_p100_beam_results.Rda")
load(file = "Results_files/p100/beam/bdgraph_sf_p100_beam_results.Rda")
load(file = "Results_files/p100/beam/huge_sf_p100_beam_results.Rda")
load(file = "Results_files/p100/beam/hubs_p100_beam_results.Rda")

load(file = "Results_files/p100/fastGHS/random_p100_fastghs_results.Rda")
load(file = "Results_files/p100/fastGHS/bdgraph_sf_p100_fastghs_results.Rda")
load(file = "Results_files/p100/fastGHS/huge_sf_p100_fastghs_results.Rda")
load(file = "Results_files/p100/fastGHS/hubs_p100_fastghs_results.Rda")

load(file = "Results_files/p200/GHSGEM/random_p200_GEM_results.Rda")
load(file = "Results_files/p200/GHSGEM/bdgraph_sf_p200_GEM_results.Rda")
load(file = "Results_files/p200/GHSGEM/huge_sf_p200_GEM_results.Rda")
load(file = "Results_files/p200/GHSGEM/hubs_p200_GEM_results.Rda")

load(file = "Results_files/p200/glasso/random_p200_glasso_results.Rda")
load(file = "Results_files/p200/glasso/bdgraph_sf_p200_glasso_results.Rda")
load(file = "Results_files/p200/glasso/huge_sf_p200_glasso_results.Rda")
load(file = "Results_files/p200/glasso/hubs_p200_glasso_results.Rda")

load(file = "Results_files/p200/beam/random_p200_beam_results.Rda")
load(file = "Results_files/p200/beam/bdgraph_sf_p200_beam_results.Rda")
load(file = "Results_files/p200/beam/huge_sf_p200_beam_results.Rda")
load(file = "Results_files/p200/beam/hubs_p200_beam_results.Rda")

scores <- c("MCC", "TPR", "FPR", "FDR", "sl_omega", "F_norm_omega_rel", "time")

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


#######
# Calculate relative F-norms.
# p = 100

GEM_random <- GEM_bdgraph <- GEM_huge <- GEM_hubs  <- c()
for (i in 1:50) {
  GEM_random[i] <- random_p100_results$results[i, "F_norm_omega"] / norm(sim_random[[i]]$omega, type = "f")
  GEM_bdgraph[i] <- bdgraph_sf_p100_results$results[i, "F_norm_omega"] / norm(sim_bdgraph_sf[[i]]$omega, type = "f")
  GEM_huge[i] <- huge_sf_p100_results$results[i, "F_norm_omega"] / norm(sim_huge_sf[[i]]$omega, type = "f")
  GEM_hubs[i] <- hubs_p100_results$results[i, "F_norm_omega"] / norm(sim_hubs[[i]]$omega, type = "f")
}

GLASSO_random <- GLASSO_bdgraph <- GLASSO_huge <- GLASSO_hubs  <- c()
for (i in 1:50) {
  GLASSO_random[i] <- random_p100_glasso_results$results[i, "F_norm_omega"] / norm(sim_random[[i]]$omega, type = "f")
  GLASSO_bdgraph[i] <- bdgraph_sf_p100_glasso_results$results[i, "F_norm_omega"] / norm(sim_bdgraph_sf[[i]]$omega, type = "f")
  GLASSO_huge[i] <- huge_sf_p100_glasso_results$results[i, "F_norm_omega"] / norm(sim_huge_sf[[i]]$omega, type = "f")
  GLASSO_hubs[i] <- hubs_p100_glasso_results$results[i, "F_norm_omega"] / norm(sim_hubs[[i]]$omega, type = "f")
}


round(cbind(rbind(mean(GEM_random), mean(GEM_bdgraph), mean(GEM_huge), mean(GEM_hubs)),
            rbind(sd(GEM_random), sd(GEM_bdgraph), sd(GEM_huge), sd(GEM_hubs))), 3)

round(cbind(rbind(mean(GLASSO_random), mean(GLASSO_bdgraph), mean(GLASSO_huge), mean(GLASSO_hubs)),
            rbind(sd(GLASSO_random), sd(GLASSO_bdgraph), sd(GLASSO_huge), sd(GLASSO_hubs))), 3)


# p = 200

GEM_p200_random <- GEM_p200_bdgraph <- GEM_p200_huge <- GEM_p200_hubs  <- c()
for (i in 1:50) {
  GEM_p200_random[i] <- random_p200_results$results[i, "F_norm_omega"] / norm(sim_random[[i]]$omega, type = "f")
  GEM_p200_bdgraph[i] <- bdgraph_sf_p200_results$results[i, "F_norm_omega"] / norm(sim_bdgraph_sf[[i]]$omega, type = "f")
  GEM_p200_huge[i] <- huge_sf_p200_results$results[i, "F_norm_omega"] / norm(sim_huge_sf[[i]]$omega, type = "f")
  GEM_p200_hubs[i] <- hubs_p200_results$results[i, "F_norm_omega"] / norm(sim_hubs[[i]]$omega, type = "f")
}

GLASSO_p200_random <- GLASSO_p200_bdgraph <- GLASSO_p200_huge <- GLASSO_p200_hubs  <- c()
for (i in 1:50) {
  GLASSO_p200_random[i] <- random_p200_glasso_results$results[i, "F_norm_omega"] / norm(sim_random[[i]]$omega, type = "f")
  GLASSO_p200_bdgraph[i] <- bdgraph_sf_p200_glasso_results$results[i, "F_norm_omega"] / norm(sim_bdgraph_sf[[i]]$omega, type = "f")
  GLASSO_p200_huge[i] <- huge_sf_p200_glasso_results$results[i, "F_norm_omega"] / norm(sim_huge_sf[[i]]$omega, type = "f")
  GLASSO_p200_hubs[i] <- hubs_p200_glasso_results$results[i, "F_norm_omega"] / norm(sim_hubs[[i]]$omega, type = "f")
}


round(cbind(rbind(mean(GEM_p200_random), mean(GEM_p200_bdgraph), mean(GEM_p200_huge), mean(GEM_p200_hubs)),
            rbind(sd(GEM_p200_random), sd(GEM_p200_bdgraph), sd(GEM_p200_huge), sd(GEM_p200_hubs))), 3)

round(cbind(rbind(mean(GLASSO_p200_random), mean(GLASSO_p200_bdgraph), mean(GLASSO_p200_huge), mean(GLASSO_p200_hubs)),
            rbind(sd(GLASSO_p200_random), sd(GLASSO_p200_bdgraph), sd(GLASSO_p200_huge), sd(GLASSO_p200_hubs))), 3)


round(colMeans(random_p100_fastghs_results$results[,scores]), 4)
round(colMeans(bdgraph_sf_p100_fastghs_results$results[,scores]), 4)
round(colMeans(huge_sf_p100_fastghs_results$results[,scores]), 4)
round(colMeans(hubs_p100_fastghs_results$results[,scores]), 4)

