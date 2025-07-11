library(GHSCM)
library(fastGHS)
library(huge)
library(doParallel)

# Set sample sizes.
sample_sizes <- c(seq(70, 220, 10), seq(250, 500, 25), seq(550, 800, 50),
                  seq(900, 1500, 100))
# Set number of variables.
p <- 150
# Set number of datasets.
n_datasets <- 20

# Run analysis if not yet done.
if(!file.exists("Results_files/fastGSH_vs_GHSCM_n70-n1500_p150.Rds")) {
  # Set inital seed for structure simulation.
  set.seed(20250313)
  # Generate network structure.
  sim <- huge.generator(n = 50, d = p , graph = "scale-free")
  # Run both fastGHS and GHS CM over sample sizes and replications using parallel computing.
  cl <- makeCluster(detectCores(logical = FALSE))
  registerDoParallel(cl)
  results <- foreach(i = 1:length(sample_sizes), .combine = "rbind", .packages = c("GHSCM", "fastGHS", "huge")) %dopar% {
    scores <- data.frame()
    for (r in 1:n_datasets) {
      n <- sample_sizes[i]
      set.seed(20250312 + n * p + i + r)
      data <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = sim$sigma)
      # Results using fastGHS.
      result <- fastGHS(data, AIC_selection = FALSE, fix_tau = TRUE, tau_sq = 5, epsilon = 1e-3)
      theta_est <- result$theta
      theta_est[abs(theta_est) < 1e-5] <- 0
      theta_est[theta_est != 0] <- 1
      diag(theta_est) <- 0
      method <- "fastGHS"
      edge_count <- sum(theta_est) / 2
      omega_est <- result$theta
      score <- calculate_scores(conf_matrix(sim$theta, theta_est))
      f_norm_rel <- norm(sim$omega - omega_est, type = "f") / norm(sim$omega, type = "f")
      sl_omega <- stein_loss(sim$omega, omega_est)
      scores <- rbind(scores, cbind(score, f_norm_rel, sl_omega, edge_count, n, r, method))
      # Results using GHS CM.
      map <- GHS_MAP_estimation(data, verbose = 0, max_iterations = 1000, tol = 1e-4)
      theta_est <- map$Theta_est
      method <- "GHSCM"
      edge_count <- sum(theta_est) / 2
      omega_est <- result$theta
      score <- calculate_scores(conf_matrix(sim$theta, theta_est))
      f_norm_rel <- norm(sim$omega - omega_est, type = "f") / norm(sim$omega, type = "f")
      sl_omega <- stein_loss(sim$omega, omega_est)
      scores <- rbind(scores, cbind(score, f_norm_rel, sl_omega, edge_count, n, r, method))
    }
    scores
  }
  stopCluster(cl)
  stopImplicitCluster()
  
  # Save results.
  saveRDS(results, "Results_files/fastGSH_vs_GHSCM_n70-n1500_p150.Rds")
}

# Load results if saved earlier.
if (file.exists("Results_files/fastGSH_vs_GHSCM_n70-n1500_p150.Rds")) {
  results <- readRDS("Results_files/fastGSH_vs_GHSCM_n70-n1500_p150.Rds")
}

# Separate fastGHS and GHS CM results into own data frames.
fastGHS_results <- results[results$method == "fastGHS", ]
GHSCM_results <- results[results$method == "GHSCM", ]

# Calculate mean over the dataset replicas.
fastGHS_agg_results <- aggregate(fastGHS_results[, 1:22], by = list(fastGHS_results$n), FUN = mean)
GHSCM_agg_results <- aggregate(GHSCM_results[, 1:22], by = list(GHSCM_results$n), FUN = mean)

# Create figures.
# Number of connections.
pdf("Figures/Supporting_information/fastGHS_and_GHSCM_n70-1500_connection_numbers.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
par(mgp = c(2, 1, 0))
par(mar = c(3.1, 3.1, 0.2, 0.2))
plot(sample_sizes, fastGHS_agg_results$edge_count, type = "l", log = "x",
     ylim = c(1, max(fastGHS_agg_results$edge_count)),
     xlab = "Sample size", ylab = "Number of connections", lty = 1)
lines(sample_sizes, fastGHS_agg_results$FPR*((p^2 - p) / 2), col = "red", lty = 4)
lines(sample_sizes, fastGHS_agg_results$TPR * (p-1), col = "blue", lty = 2)
abline(h = p-1, lty = 3, lwd = 2, col = "black")
legend("topright", c("All estimated", "False positives", "True positives",
                     "In true network"),
       col = c("black", "red", "blue", "black"), lty = c(1,4,2,3),
       lwd = c(1,1,1,2))
grid()

plot(sample_sizes, GHSCM_agg_results$edge_count, type = "l", log = "x",
     ylim = c(0, max(GHSCM_agg_results$edge_count)),
     xlab = "Sample size", ylab = "Number of connections", lty = 1)
lines(sample_sizes, GHSCM_agg_results$FPR*((p^2 - p) / 2), col = "red", lty = 4)
lines(sample_sizes, GHSCM_agg_results$TPR * (p-1), col = "blue", lty = 2)
abline(h = p-1, lty = 3, lwd = 2, col = "black")
grid()
dev.off()

# MCCs and F1-scores.
pdf("Figures/Supporting_information/fastGHS_and_GHSCM_n70-1500_MCC_and_F1.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
par(mgp = c(2, 1, 0))
par(mar = c(3.1, 3.1, 0.2, 0.2))
plot(sample_sizes, fastGHS_agg_results$MCC, type = "l", log = "x",
     ylim = c(0, 1), xlab = "Sample size", ylab = "MCC")
lines(sample_sizes, GHSCM_agg_results$MCC, lty = 2, lwd = 2)
legend("bottomright", legend = c("fastGHS", "GHS CM"), lty = c(1,2), lwd = c(1,2))
grid()

plot(sample_sizes, fastGHS_agg_results$F1, type = "l", log = "x", ylim = c(0, 1),
     xlab = "Sample size", ylab = expression(paste(F[1], "-score", sep = "")))
lines(sample_sizes, GHSCM_agg_results$F1, lty = 2, lwd = 2)
grid()
dev.off()

# Print F1-scores and MCCs.
data.frame(n = sample_sizes, fastGHS = fastGHS_agg_results$F1, GHSCM = GHSCM_agg_results$F1)
data.frame(n = sample_sizes, fastGHS = fastGHS_agg_results$MCC, GHSCM = GHSCM_agg_results$MCC)


