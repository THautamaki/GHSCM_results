run_analysis <- function(structure, n, p, method = "GHSGEM", p0 = 0, save_results = TRUE) {
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
  if (tolower(method) == "ghsgem") {
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
    stop("Wrong method! Possible choices are 'GHSGEM', 'GLASSO', 'pulsar' and 'fastGHS'.")
  }
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  start <- Sys.time()
  results <- foreach (i = 1:n_datasets, .combine = "rbind", .packages = packages,
                      .verbose = FALSE) %dopar% {
    sim <- simulations[[i]]
    if (tolower(method) == "ghsgem") {
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
      lmax <- getMaxCov(cov(sim$data))
      lams <- getLamPath(lmax, lmax * 0.1, len = 50)
      #lams <- huge(sim$data, nlambda = 50, method = "glasso")$lambda
      farguments <- list(lambda = lams, method = "glasso", verbose = FALSE)
      out.p <- pulsar(sim$data, fun = huge, fargs = farguments, rep.num = 20, criterion = 'stars',
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
    f_norm_sigma <- NA
    f_norm_omega <- NA
    f_norm_rel <- NA
    if (tolower(method) == "ghsgem" | tolower(method) == "fastghs") {
      sl_omega <- stein_loss(sim$omega, omega_est)
      sl_sigma <- stein_loss(sim$sigma, sigma_est)
      f_norm_omega <- norm(sim$omega - omega_est, type = "f")
      f_norm_sigma <- norm(sim$sigma - sigma_est, type = "f")
      f_norm_rel <- f_norm_omega / norm(sim$omega, type = "f")
    }
    else if (tolower(method) == "glasso") {
      sl_omega <- stein_loss(sim$omega, omega_est)
      f_norm_omega <- norm(sim$omega - omega_est, type = "f")
      f_norm_rel <- f_norm_omega / norm(sim$omega, type = "f")
    }
    score <- cbind(i, calculate_scores(cm), edge_count, sl_sigma, sl_omega, f_norm_sigma,
                   f_norm_omega, f_norm_rel, time)
    score
  }
  total_time <- Sys.time() - start
  stopCluster(cl)
  stopImplicitCluster()
  results <- list(results = results, total_time = total_time)
  if (save_results) {
    results_path <- paste0("Results_files/p", p, "/", method, "/")
    if (tolower(method) == "ghsgem" & tolower(structure) == "random" & p0 > 0) {
      saveRDS(results, file = paste0(results_path, structure, "_p", p, "_", method, "_results_2.Rds"))
    }
    else {
      saveRDS(results, file = paste0(results_path, structure, "_p", p, "_", method, "_results.Rds"))
    }
  }
  return(list(results = results, total_time = total_time))
}

run_multiple_methods <- function(methods, structures, sample_sizes, variable_numbers, save_results = TRUE) {
  all_results <- list()
  for (method in methods) {
    for (structure in structures) {
      for (n in sample_sizes) {
        for (p in variable_numbers) {
          all_results[[method]][[structure]][[paste0("n", n, "_p", p)]] <- run_analysis(structure, n, p, method, save_results = save_results)
          if (method == "GHSGEM" & structure == "random") {
            all_results[[paste0(method, "_p/2")]][[structure]][[paste0("n", n, "_p", p)]] <- run_analysis(structure, n, p, method, p0 = p/2, save_results = save_results)
          }
        }
      }
    }
  }
  return(all_results)
}

add_R_results <- function(results_list, methods, structures, sample_sizes, variable_numbers) {
  for (method in methods) {
    for (structure in structures) {
      for (n in sample_sizes) {
        for (p in variable_numbers) {
          path <- paste0("Results_files/p",p, "/", method, "/")
          if (tolower(method) == "ghsgem" & tolower(structure) == "random") {
            results <- readRDS(paste0(path, structure, "_p", p, "_", method, "_results_2.Rds"))
            results_list[[paste0(method, "_p/2")]][[structure]][[paste0("n", n, "_p", p)]]$results <- results$results
            results_list[[paste0(method, "_p/2")]][[structure]][[paste0("n", n, "_p", p)]]$total_time <- results$total_time
          }
          path <- paste0("Results_files/p",p, "/", method, "/")
          results <- readRDS(paste0(path, structure, "_p", p, "_", method, "_results.Rds"))
          results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results <- results$results
          results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$total_time <- results$total_time
        }
      }
    }
  }
  return(results_list)
}

add_MATLAB_results <- function(results_list, methods, structures, sample_sizes, variable_numbers) {
  for (method in methods) {
    for (structure in structures) {
      for (n in sample_sizes) {
        for (p in variable_numbers) {
          path <- paste0("Results_files/p",p, "/", method, "/")
          results <- read.csv(paste0(path, structure, "_p", p, "_scores.csv"))
          times <- read.csv(paste0("Results_files/", method, "_total_times.csv"))
          time <- times[times$p == p & times$structure == structure, "time"]
          results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results <- results
          results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$total_time <- time
          if (method == "GHS_LLA") {
            alg_time <- results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results$time
            tau_time <- results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results$tau_time
            tau_prop <- tau_time / (alg_time + tau_time)
            results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results$time <- alg_time + tau_time
            results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results$alg_time <- alg_time
            results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results$tau_prop <- tau_prop
          }
        }
      }
    }
  }
  return(results_list)
}

print_results <- function(scores, results_list, methods, structure, n, p) {
  mean_results <- sd_results <- data.frame()
  for (method in methods) {
    mean_scores <- apply(results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results[, scores],
                         2, mean)
    mean_scores <- as.list(mean_scores)
    tot_time <- results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$total_time
    mean_scores$total_time <- as.double(tot_time, units = "secs")
    sd_scores <- apply(results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results[, scores],
                       2, sd)
    mean_results <- rbind(mean_results, mean_scores)
    sd_results <- rbind(sd_results, sd_scores)
  }
  rownames(mean_results) <- rownames(sd_results) <- methods
  colnames(sd_results) <- scores
  cat(paste0("Mean of the performance scores with the ", structure, " network structure\n"))
  print(round(mean_results, 4))
  cat("\nSd of the performance scores\n")
  print(round(sd_results, 4))
}

create_latex_table <- function(scores, results_list, methods, structure, n, p) {
  options(scipen = 999)
  for (method in methods) {
    cat(format(method, width = 14), " & ", sep = "")
    for (score in scores) {
      score_mean <- mean(results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results[, score])
      score_sd <- sd(results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results[, score])
      if (score == "FPR") {
        cat(format(paste0(sprintf("%.4f", round(score_mean, 4)), " (",
                          sprintf("%.4f", round(score_sd, 4)), ") & "),  width = 24))
      }
      else if (score == "time") {
        tot_time <- results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$total_time
        tot_time <- as.double(tot_time, units = "secs")
        rf_tot <- select_rounding(tot_time)
        rf_mean <- select_rounding(score_mean)
        cat(paste0(format(sprintf(rf_mean[2], round(score_mean, as.numeric(rf_mean[1]))), width = 12),
                   " & ", sprintf(rf_tot[2], round(tot_time, as.numeric(rf_tot[1]))), " \\\\ \n"))
      }
      else {
        if (is.na(score_mean)) {
          cat(format(".. & ", width = 22))
        }
        else {
          cat(format(paste0(sprintf("%.3f", round(score_mean, 3)), " (",
                            sprintf("%.3f", round(score_sd, 3)), ") & "), width = 22))
        }
      }
    }
  }
  options(scipen = 0)
}

select_rounding <- function(value) {
  if (value < 1) {
    rounding <- 3
    format <- "%.3f"
  }
  else if (value < 10) {
    rounding <- 2
    format <- "%.2f"
  }
  else if (value < 100) {
    rounding <- 1
    format <- "%.1f"
  }
  else {
    rounding <- 0
    format <- "%.0f"
  }
  return(c(rounding, format))
}
