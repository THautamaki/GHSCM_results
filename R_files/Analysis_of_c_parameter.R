library(GHSGEM)
library(huge)
library(doParallel)

# Set number of variables, sample sizes, number of datasets and target FDR.
p <- c(seq(100, 250, 10), seq(275, 350, 25), c(400, 500))
n <- p - 10
n_repeats <- 50
target_fdr <- 0.2

# Run analysis with FDR-controlled c if not yet done.
if(!file.exists("Results_files/Parameter_c/Results_with_fdr_control.Rds")) {
  cl <- makeCluster(detectCores(logical = FALSE))
  registerDoParallel(cl)
  (start <- Sys.time())
  results_with_fdr_control <- foreach (r = 1:n_repeats, .combine = "rbind",
                                       .packages = c("GHSGEM", "huge"), .verbose = TRUE) %dopar% {
    scores <- data.frame()
    for (i in 1:length(p)) {
      n1 <- n[i]
      p1 <- p[i]
      p0 <- p1 - 1
      if (p1 < 250) {
        c <- 50
        c_vec <- 1:100
      }
      else {
        c <- p1/5
        c_vec <- (c-50):(c+50)
      }
      tau_f <- (p0/(p1 * (p1 - 1)/2)) * (c * sqrt(p0)/(n1 * sqrt(n1)))
      set.seed(20250322 + p1 * n1 + r)
      sim <- huge.generator(n = n1, d = p1, graph = "scale-free")
      map <- graphical_horseshoe_map_Cpp(sim$data, fixed_tau = tau_f, max_iter = 500)
      theta_est <- map$Theta_est
      cm <- conf_matrix(sim$theta, theta_est)
      fdr <- calculate_scores(cm)$FDR
      if (is.na(fdr)) {
        fdr <- 0
      }
      iter <- 1
      # Use "brute force" for smaller problems.
      if (p1 < 170) {
        while (abs(target_fdr - fdr) > 0.02 & iter < 50) {
          if (target_fdr - fdr < 0) {
            c <- c - 1
          }
          else {
            c <- c + 1
          }
          tau_f <- (p0/(p1 * (p1 - 1)/2)) * (c * sqrt(p0)/(n1 * sqrt(n1)))
          map <- graphical_horseshoe_map_Cpp(sim$data, fixed_tau = tau_f, max_iter = 500)
          theta_est <- map$Theta_est
          cm <- conf_matrix(sim$theta, theta_est)
          fdr <- calculate_scores(cm)$FDR
          if (is.na(fdr)) {
            fdr <- 0
          }
          iter <- iter + 1
        }
      }
      # Use binary search otherwise.
      else {
        while (abs(target_fdr - fdr) > 0.02 & iter < 8) {
          if (target_fdr - fdr < 0) {
            c_vec <- c_vec[1:floor(length(c_vec)/2)]
            c <- c_vec[floor(length(c_vec)/2)]
          }
          else {
            c_vec <- c_vec[floor(length(c_vec)/2):length(c_vec)]
            c <- c_vec[floor(length(c_vec)/2)]
          }
          tau_f <- (p0/(p1 * (p1 - 1)/2)) * (c * sqrt(p0)/(n1 * sqrt(n1)))
          map <- graphical_horseshoe_map_Cpp(sim$data, fixed_tau = tau_f, max_iter = 500)
          theta_est <- map$Theta_est
          cm <- conf_matrix(sim$theta, theta_est)
          fdr <- calculate_scores(cm)$FDR
          if (is.na(fdr)) {
            fdr <- 0
          }
          iter <- iter + 1
        }
      }
      score <- cbind(p1, n1, r, calculate_scores(cm), tau_f, c, fdr, iter)
      scores <- rbind(scores, score)
    }
    scores
  }
  end <- Sys.time()
  (end - start)
  stopCluster(cl)
  stopImplicitCluster()
  # Save results.
  saveRDS(results_with_fdr_control, "Results_files/Parameter_c/Results_with_fdr_control.Rds")
}

# Run analysis with default c if not yet done.
if(!file.exists("Results_files/Parameter_c/Results_with_defaults.Rds")) {
  cl <- makeCluster(detectCores(logical = FALSE))
  registerDoParallel(cl)
  (start <- Sys.time())
  results_with_defaults <- foreach (r = 1:n_repeats, .combine = "rbind",
                                    .packages = c("GHSGEM", "huge"), .verbose = TRUE) %dopar% {
    scores <- data.frame()
    for (i in 1:length(p)) {
      n1 <- n[i]
      p1 <- p[i]
      set.seed(20250322 + p1 * n1 + r)
      sim <- huge.generator(n = n1, d = p1, graph = "scale-free")
      map <- GHS_MAP_estimation(sim$data, fixed_tau = tau_f, max_iter = 500)
      theta_est <- map$Theta_est
      cm <- conf_matrix(sim$theta, theta_est)
      score <- cbind(p1, n1, r, calculate_scores(cm))
      scores <- rbind(scores, score)
    }
    scores
  }
  end <- Sys.time()
  (end - start)
  stopCluster(cl)
  stopImplicitCluster()
  # Save results.
  saveRDS(results_with_defaults, "Results_files/Parameter_c/Results_with_defaults.Rds")
}

# Load FDR-controlled results if saved earlier.
if (file.exists("Results_files/Parameter_c/Results_with_fdr_control.Rds")) {
  results_with_fdr_control <- readRDS("Results_files/Parameter_c/Results_with_fdr_control.Rds")
}

# Load results with defaults if saved earlier.
if (file.exists("Results_files/Parameter_c/Results_with_defaults.Rds")) {
  results_with_defaults <- readRDS("Results_files/Parameter_c/Results_with_defaults.Rds")
}

# Add number of edges into results.
results_with_fdr_control$edge_count <- (results_with_fdr_control$TPR + results_with_fdr_control$FPR) * 
  (results_with_fdr_control$p1 - 1)
results_with_defaults$edge_count <- (results_with_defaults$TPR + results_with_defaults$FPR) * 
  (results_with_defaults$p1 - 1)

# Aggregate results over p.
agg_results <- aggregate(results_with_fdr_control, by = list(results_with_fdr_control$p1),
                         FUN = mean, na.rm = TRUE)
agg_results2 <- aggregate(results_with_defaults, by = list(results_with_defaults$p1), FUN = mean)

# Calculate what denominator would be instead of 5. Print denominator, mean of c and difference of
# mean of c and p / 5.
data.frame(p = agg_results$p1, denominator = agg_results$p1 / agg_results$c, c_orig = agg_results$p1/5, 
           c = agg_results$c, difference = agg_results$c - agg_results$p1/5)

# Plot key performance scores.
pdf("Figures/Supplementary/c_parameter_scores.pdf", width = 7, height = 9)
par(mfrow = c(3,2))
par(mgp = c(2, 1, 0))
par(mar = c(3.1, 3.1, 0.2, 0.2))
lty <- 2
col = "blue"
score <- "FDR"
plot(p, agg_results[,score], type = "l", xlab = "Number of variables", ylab = score,
     ylim = c(min(agg_results[,score], agg_results[,score]),
              max(agg_results[,score], agg_results2[,score])),
     col = col, lty = lty)
lines(p, agg_results2[,score], col = "black", lty = 1, lwd = 1)
legend("topright", legend = c("Default c", "FDR-controlled c"), lty = c(1, lty), col = c("black", col),
       lwd = c(1,1))
grid()

score <- "FPR"
plot(p, agg_results[,score], type = "l", xlab = "Number of variables", ylab = score,
     ylim = c(0, max(agg_results[,score], agg_results2[,score])),
     col = col, lty = lty)
lines(p, agg_results2[,score], col = "black", lty = 1, lwd = 1)
grid()

score <- "MCC"
plot(p, agg_results[,score], type = "l", xlab = "Number of variables", ylab = score,
     ylim = c(min(agg_results[,score], agg_results[,score]),
              max(agg_results[,score], agg_results2[,score])),
     col = col, lty = lty)
lines(p, agg_results2[,score], col = "black", lty = 1, lwd = 1)
grid()

score <- "F1"
plot(p, agg_results[,score], type = "l", ylim = c(min(agg_results[,score], agg_results[,score]),
                                                  max(agg_results[,score], agg_results2[,score])),
     xlab = "Number of variables", ylab = score, col = col, lty = lty)
lines(p, agg_results2[,score], col = "black", lty = 1, lwd = 1)
grid()

score <- "TPR"
plot(p, agg_results[,score], type = "l", ylim = c(min(agg_results[,score], agg_results[,score]),
                                                  max(agg_results[,score], agg_results2[,score])),
     xlab = "Number of variables", ylab = score, col = col, lty = lty)
lines(p, agg_results2[,score], col = "black", lty = 1, lwd = 1)
grid()

score <- "Total number of connections"
plot(p, agg_results$edge_count, type = "l",
     ylim = c(0, max(agg_results$edge_count, agg_results2$edge_count)),
     xlab = "Number of variables", ylab = score, col = col, lty = lty)
lines(p, agg_results2$edge_count, col = "black", lty = 1, lwd = 1)
grid()

dev.off()
