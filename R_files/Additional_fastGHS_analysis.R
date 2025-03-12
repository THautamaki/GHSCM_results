library(GHSGEM)
library(fastGHS)
library(huge)
library(doParallel)

sample_sizes <- c(seq(70, 220, 10), seq(250, 500, 25))
p <- 150

cl <- makeCluster(detect_cores(logical = FALSE))
registerDoParallel(cl)
results <- foreach(i = 1:length(sample_sizes), .combine = "rbind", .packages = c("GHSGEM", "fastGHS", "huge")) %dopar% {
  scores <- data.frame()
  for (r in 1:10) {
    n <- sample_sizes[i]
    set.seed(20250312 + n * p + i + r)
    sim <- huge.generator(n = n, d = p , graph = "scale-free")
    # Results using fastGHS.
    result <- fastGHS(sim$data, AIC_selection = FALSE, fix_tau = TRUE, tau_sq = 5, epsilon = 1e-3)
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
    # Results using GHS GEM.
    map <- GHS_MAP_estimation(sim$data, verbose = 0, max_iterations = 1000, tol = 1e-4)
    theta_est <- map$Theta_est
    method <- "GHSGEM"
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

# Separate fastGHS and GHS GEM results into own data frames.
fastGHS_results <- results[results$method == "fastGHS", ]
GHSGEM_results <- results[results$method == "GHSGEM", ]

# Calculate mean over the dataset replicas.
fastGHS_agg_results <- aggregate(fastGHS_results[, 1:22], by = list(fastGHS_results$n), FUN = mean)
GHSGEM_agg_results <- aggregate(GHSGEM_results[, 1:22], by = list(GHSGEM_results$n), FUN = mean)


pdf("Figures/Supplementary/fastGHS_different_sample_sizes.pdf", width = 9, height = 8.1)
par(mar = c(4.1, 4.1, 0.1, 0.1))
plot(sample_sizes, fastGHS_agg_results$edge_count, type = "l", log = "x",
     ylim = c(0, max(fastGHS_agg_results$edge_count)),
     xlab = "Number of observations", ylab = "Number of connections")
lines(sample_sizes, fastGHS_agg_results$FPR*((p^2 - p) / 2), col = "blue")
lines(sample_sizes, fastGHS_agg_results$TPR * (p-1), col = "green")

lines(sample_sizes, GHSGEM_agg_results$edge_count, lwd = 2, lty = 2)
lines(sample_sizes, GHSGEM_agg_results$FPR*((p^2 - p) / 2), col = "blue", lwd = 2, lty = 2)
lines(sample_sizes, GHSGEM_agg_results$TPR * (p-1), col = "green", lwd = 2, lty = 2)

abline(h = p-1, lty = 2, col = "black")
legend("topright", c("fastGHS: all connections", "fastGHS: false connections", "fastGHS: correct connections",
                     "GHS GEM: all connections", "GHS GEM: false connections", "GHS GEM: correct connections",
                     "Number of true connections"),
       col = c("black", "blue", "green", "black", "blue", "green", "black"), lty = c(1,1,1,2,2,2,2),
       lwd = c(1,1,1,2,2,2,1))
grid()
dev.off()
