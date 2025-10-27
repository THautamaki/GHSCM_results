run_analysis <- function(structure, n, p, method = "GHSCM", p0 = 0, save_results = TRUE) {
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
  if (tolower(method) == "ghscm") {
    packages <- c("GHSCM")
  }
  else if (tolower(method) == "glasso") {
    packages <- c("huge", "GHSCM")
  }
  else if (tolower(method) == "beam") {
    packages <- c("beam", "GHSCM")
  }
  else if (tolower(method) == "pulsar") {
    packages <- c("pulsar", "huge", "GHSCM")
  }
  else if (tolower(method) == "fastghs") {
    packages <- c("fastGHS", "GHSCM")
  }
  else {
    stop("Wrong method! Possible choices are 'GHSCM', 'GLASSO', 'pulsar', 'beam', and 'fastGHS'.")
  }
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  start <- Sys.time()
  results <- foreach (i = 1:n_datasets, .combine = "rbind", .packages = packages,
                      .verbose = FALSE) %dopar% {
    sim <- simulations[[i]]
    if (tolower(method) == "ghscm") {
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
    if (tolower(method) == "ghscm" | tolower(method) == "fastghs") {
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
    if (tolower(method) == "ghscm" & tolower(structure) == "random" & p0 > 0) {
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
          if (method == "GHSCM" & structure == "random") {
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
          if (tolower(method) == "ghscm" & tolower(structure) == "random") {
            results <- readRDS(paste0(path, structure, "_p", p, "_", method, "_results_2.Rds"))
            results_list[[paste0(method, "_p/2")]][[structure]][[paste0("n", n, "_p", p)]]$results <- results$results
            results_list[[paste0(method, "_p/2")]][[structure]][[paste0("n", n, "_p", p)]]$total_time <- as.double(results$total_time, units = "secs")
          }
          results <- readRDS(paste0(path, structure, "_p", p, "_", method, "_results.Rds"))
          results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results <- results$results
          results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$total_time <- as.double(results$total_time, units = "secs")
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

print_results <- function(scores, results_list, methods, structure, n, p, return = FALSE) {
  mean_results <- sd_results <- data.frame()
  add_total_time <- FALSE
  if (any(scores == "total_time")) {
    scores <- scores[-which(scores == "total_time")]
    add_total_time <- TRUE
  }
  for (method in methods) {
    if (length(scores) > 1) {
      mean_scores <- apply(results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results[, scores],
                           2, mean)
      sd_scores <- apply(results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results[, scores],
                         2, sd)
    }
    else {
      mean_scores <- apply(as.data.frame(results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results[, scores]),
                           2, mean)
      sd_scores <- apply(as.data.frame(results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results[, scores]),
                         2, sd)
    }
    names(mean_scores) <- names(sd_scores) <- scores
    mean_scores <- as.list(mean_scores)
    sd_scores <- as.list(sd_scores)
    if (add_total_time) {
      tot_time <- results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$total_time
      mean_scores$total_time <- as.double(tot_time, units = "secs")
      sd_scores$total_time <- NA
    }
    mean_results <- rbind(mean_results, mean_scores)
    sd_results <- rbind(sd_results, sd_scores)
  }
  rownames(mean_results) <- rownames(sd_results) <- methods
  if (return) {
    return(list(means = mean_results, sds = sd_results))
  }
  cat(paste0("Mean of the performance scores with the ", structure, " network structure\n"))
  print(round(mean_results, 4))
  cat("\nSd of the performance scores\n")
  print(round(sd_results, 4))
}

create_latex_table <- function(scores, results_list, methods, structure, n, p, highlight_best = TRUE) {
  options(scipen = 999)
  means_and_sds <- print_results(scores, results_list, methods, structure, n, p, TRUE)
  df_means <- means_and_sds$means
  df_sds <- means_and_sds$sds
  for (method in methods) {
    if (method == "GHSCM") {
      cat(format("GHS CM, $p_0=p-1$", width = 19), sep = "")
    }
    else if (method == "GHSCM_p/2") {
      cat(format("GHS CM, $p_0=p/2$", width = 19), sep = "")
    }
    else if (method == "GHS_MCMC") {
      cat(format("GHS MCMC", width = 19), sep = "")
    }
    else if (method == "GHS_LLA") {
      cat(format("GHS LLA", width = 19), sep = "")
    }
    else if (method == "HSL_MCMC") {
      cat(format("GHS-like MCMC", width = 19), sep = "")
    }
    else if (method == "HSL_ECM") {
      cat(format("GHS-like ECM", width = 19), sep = "")
    }
    else if (method == "GLASSO") {
      cat(format("GLASSO (StARS)", width = 19), sep = "")
    }
    else if (method == "fastGHS") {
      cat(format("GHS ECM", width = 19), sep = "")
    }
    else {
      cat(format(method, width = 19), sep = "")
    }
    for (score in scores) {
      #score_mean <- mean(results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results[, score])
      #score_sd <- sd(results_list[[method]][[structure]][[paste0("n", n, "_p", p)]]$results[, score])
      score_mean <- df_means[method, score]
      score_sd <- df_sds[method, score]
      best <- FALSE
      if (!is.na(score_mean)) {
        decr <- TRUE
        if (any(score == c("FPR", "FDR", "sl_omega", "f_norm_rel"))) {
          decr <- FALSE
        }
        ordered_scores <- df_means[, score][order(df_means[, score], decreasing = decr)]
        if (any(score_mean == ordered_scores[1:2])) {
          best <- TRUE
        }
      }
      if (score == "FPR") {
        mean_score <- sprintf("%.4f", round(score_mean, 4))
        if (best & highlight_best) {
          mean_score <- paste0("\\textbf{", sprintf("%.4f", round(score_mean, 4)), "}")
        }
        cat(format(paste0(" & ", mean_score, " (",
                          sprintf("%.4f", round(score_sd, 4)), ")"),  width = 27))
      }
      else if (score == "time" | score == "total_time") {
        rf_mean <- select_rounding(score_mean)
        cat(paste0(" & ", format(sprintf(rf_mean[2], round(score_mean, as.numeric(rf_mean[1]))), width = 6)))
      }
      else if (score == "sl_omega") {
        if (is.na(score_mean)) {
          cat(format(" & .. ", width = 23))
        }
        else {
          mean_score <- sprintf("%.2f", round(score_mean, 2))
          if (best & highlight_best) {
            mean_score <- paste0("\\textbf{", sprintf("%.2f", round(score_mean, 2)), "}")
          }
          cat(format(paste0(" & ", mean_score, " (",
                            sprintf("%.2f", round(score_sd, 2)), ")"), width = 23))
        }
      }
      else {
        if (is.na(score_mean)) {
          cat(format(" & .. ", width = 25))
        }
        else {
          mean_score <- sprintf("%.3f", round(score_mean, 3))
          if (best & highlight_best) {
            mean_score <- paste0("\\textbf{", sprintf("%.3f", round(score_mean, 3)), "}")
          }
          cat(format(paste0(" & ", mean_score, " (",
                            sprintf("%.3f", round(score_sd, 3)), ")"), width = 25))
        }
      }
    }
    cat(" \\\\ \n")
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

calculate_runtimes <- function(all_results, methods, structures, sample_sizes, variable_numbers) {
  # Runtimes over all datasets and structures.
  runtimes <- data.frame()
  for (method in other_methods) {
    for (n in sample_sizes) {
      c_runtime <- c_totaltime <- c()
      for (p in variable_numbers) {
        sum_runtime <- 0
        sum_totaltime <- 0
        for (structure in structures) {
          sum_runtime <- sum_runtime + sum(all_results[[method]][[structure]][[paste0("n", n, "_p", p)]]$results$time)
          sum_totaltime <- sum_totaltime + all_results[[method]][[structure]][[paste0("n", n, "_p", p)]]$total_time
        }
        c_runtime <- c(c_runtime, sum_runtime)
        c_totaltime <- c(c_totaltime, sum_totaltime)
      }
    }
    runtimes <- rbind(runtimes, c(c_runtime / 200, c_totaltime / 4))
  }
  colnames(runtimes) <- c("rt_p100", "rt_p200", "tt_p100", "tt_p200")
  rownames(runtimes) <- other_methods
  return(runtimes)
}

create_runtime_table <- function(all_results, methods, structures, sample_sizes, variable_numbers) {
  runtimes <- calculate_runtimes(all_results, methods, structures, sample_sizes, variable_numbers)
  for (method in methods) {
    if (method == "GHSCM") {
      cat(format("GHS CM, $p_0=p-1$", width = 19), sep = "")
    }
    else if (method == "GHSCM_p/2") {
      cat(format("GHS CM, $p_0=p/2$", width = 19), sep = "")
    }
    else if (method == "GHS_MCMC") {
      cat(format("GHS MCMC", width = 19), sep = "")
    }
    else if (method == "GHS_LLA") {
      cat(format("GHS LLA", width = 19), sep = "")
    }
    else if (method == "HSL_MCMC") {
      cat(format("GHS-like MCMC", width = 19), sep = "")
    }
    else if (method == "HSL_ECM") {
      cat(format("GHS-like ECM", width = 19), sep = "")
    }
    else if (method == "GLASSO") {
      cat(format("GLASSO (StARS)", width = 19), sep = "")
    }
    else if (method == "fastGHS") {
      cat(format("GHS ECM", width = 19), sep = "")
    }
    else {
      cat(format(method, width = 19), sep = "")
    }
    for (time in runtimes[method,]) {
      r_and_f <- select_rounding(time)
      cat(" & ", format(sprintf(r_and_f[2], round(time, as.numeric(r_and_f[1]))), width = 12))
    }
    cat("\\\\ \n")
  }
}

calculate_false_positives <- function(all_results, methods, structures, sample_sizes, variable_numbers) {
  other_methods <- methods[-1]
  df_false_positives <- data.frame(row.names = methods)
  i <- 1
  for (n in sample_sizes) {
    for (p in variable_numbers) {
      for (structure in structures) {
        if (structure == "random") {
          positives <- print_results(c("TPR", "FPR"), all_results, methods, structure, n, p, return = TRUE)
          df_false_positives[, i] <- positives$means$FPR * (p*(p-1)/2 - (p/2))
        }
        else if (structure == "hubs") {
          positives <- print_results(c("TPR", "FPR"), all_results, other_methods, structure, n, p, return = TRUE)
          df_false_positives[other_methods, i] <- positives$means$FPR * (p*(p-1)/2 - (p - p/20))
        }
        else {
          positives <- print_results(c("TPR", "FPR"), all_results, other_methods, structure, n, p, return = TRUE)
          df_false_positives[other_methods, i] <- positives$means$FPR * (p*(p-1)/2 - (p - 1))
        }
        i <- i + 1
      }
    }
  }
  colnames(df_false_positives) <- c(structures, structures)
  return(df_false_positives)
}

create_false_positives_table <- function(all_results, methods, structures, sample_sizes, variable_numbers) {
  false_positives <- calculate_false_positives(all_results, methods, structures, sample_sizes, variable_numbers)
  for (method in methods) {
    if (method == "GHSCM") {
      cat(format("GHS CM, $p_0=p-1$", width = 19), sep = "")
    }
    else if (method == "GHSCM_p/2") {
      cat(format("GHS CM, $p_0=p/2$", width = 19), sep = "")
    }
    else if (method == "GHS_MCMC") {
      cat(format("GHS MCMC", width = 19), sep = "")
    }
    else if (method == "GHS_LLA") {
      cat(format("GHS LLA", width = 19), sep = "")
    }
    else if (method == "HSL_MCMC") {
      cat(format("GHS-like MCMC", width = 19), sep = "")
    }
    else if (method == "HSL_ECM") {
      cat(format("GHS-like ECM", width = 19), sep = "")
    }
    else if (method == "GLASSO") {
      cat(format("GLASSO (StARS)", width = 19), sep = "")
    }
    else if (method == "fastGHS") {
      cat(format("GHS ECM", width = 19), sep = "")
    }
    else {
      cat(format(method, width = 19), sep = "")
    }
    for (fp in false_positives[method,]) {
      if (is.na(fp)) {
        cat(format(" & .. ", width = 12))
        next
      }
      cat(" & ", format(sprintf("%i", round(fp, 0)), width = 9), sep = "")
    }
    cat(" \\\\ \n")
  }
}

create_sparsity_table <- function(all_results, methods, structures, sample_sizes, variable_numbers) {
  for (n in sample_sizes) {
    for (p in variable_numbers) {
      cat("p: ", p, "\n")
      for (method in random_methods) {
        if (method == "GHSCM") {
          cat(format("GHS CM, $p_0=p-1$", width = 19), sep = "")
        }
        else if (method == "GHSCM_p/2") {
          cat(format("GHS CM, $p_0=p/2$", width = 19), sep = "")
        }
        else if (method == "GHS_MCMC") {
          cat(format("GHS MCMC", width = 19), sep = "")
        }
        else if (method == "GHS_LLA") {
          cat(format("GHS LLA", width = 19), sep = "")
        }
        else if (method == "HSL_MCMC") {
          cat(format("GHS-like MCMC", width = 19), sep = "")
        }
        else if (method == "HSL_ECM") {
          cat(format("GHS-like ECM", width = 19), sep = "")
        }
        else if (method == "GLASSO") {
          cat(format("GLASSO (StARS)", width = 19), sep = "")
        }
        else if (method == "fastGHS") {
          cat(format("GHS ECM", width = 19), sep = "")
        }
        else {
          cat(format(method, width = 19), sep = "")
        }
        for (structure in structures) {
          edge_count <- all_results[[method]][[structure]][[paste0("n", n, "_p", p)]]$results$edge_count
          sparsity <- 1 - edge_count / (p * (p - 1) / 2)
          mean_sparsity <- mean(sparsity) * 100
          sd_sparsity <- sd(sparsity) * 100
          cat(" & ", sprintf("%.1f", round(mean_sparsity, 1)), " (", sprintf("%.1f", round(sd_sparsity, 1)), ") ", sep = "")
        }
        cat("\\\\\n")
      }
    }
  }
}
