library(GHSCM)

# Set data dimensions.
n <- 120
p <- 100

# Set path (no needed to change).
path <- paste0("Data/n", n, "_p", p, "/")

# Load datasets.
sim_bdgraph_sf <- readRDS(file = paste0(path, "bdgraph_scale-free_n", n, "_p", p, ".Rds"))

# Create plot for normal case example.
data_nro <- 2

map <- GHS_MAP_estimation(sim_bdgraph_sf[[data_nro]]$data, verbose = 1, max_iterations = 500)

kappa <- 1 - map$Kappa

thresholds <- seq(0.25, 0.75, 0.01)
threshold_scores <- data.frame()
for(thr in thresholds) {
  theta_est <- kappa
  theta_est[theta_est < thr] <- 0
  theta_est[theta_est != 0] <- 1
  diag(theta_est) <- 0
  scores <- calculate_scores(conf_matrix(sim_bdgraph_sf[[data_nro]]$theta, theta_est))
  threshold_scores <- rbind(threshold_scores, cbind(thr, scores))
}

pdf("Figures/Supplementary/Thresholding_normal_case_example.pdf", width = 12, height = 3)
par(mfrow = c(1,4))
par(mgp = c(2, 1, 0))
par(mar = c(3.1, 3.1, 0.2, 0.2))
plot(threshold_scores$thr, threshold_scores$MCC, xlab = "Threshold of 1 - kappa", ylab = "MCC", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "MCC")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "MCC"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$F1, xlab = "Threshold of 1 - kappa", ylab = "F1", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "F1")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "F1"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$TPR, xlab = "Threshold of 1 - kappa", ylab = "TPR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "TPR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "TPR"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$FDR, xlab = "Threshold of 1 - kappa", ylab = "FDR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "FDR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "FDR"], lty = 2, col = "blue")
grid()
dev.off()

# Create plot for worse case example.
data_nro <- 5

map <- GHS_MAP_estimation(sim_bdgraph_sf[[data_nro]]$data, verbose = 1, max_iterations = 500)

kappa <- 1 - map$Kappa

thresholds <- seq(0.25, 0.75, 0.01)
threshold_scores <- data.frame()
for(thr in thresholds) {
  theta_est <- kappa
  theta_est[theta_est < thr] <- 0
  theta_est[theta_est != 0] <- 1
  diag(theta_est) <- 0
  scores <- calculate_scores(conf_matrix(sim_bdgraph_sf[[data_nro]]$theta, theta_est))
  threshold_scores <- rbind(threshold_scores, cbind(thr, scores))
}

pdf("Figures/Supplementary/Thresholding_worse_case_example.pdf", width = 12, height = 3)
par(mfrow = c(1,4))
par(mgp = c(2, 1, 0))
par(mar = c(3.1, 3.1, 0.2, 0.2))
plot(threshold_scores$thr, threshold_scores$MCC, xlab = "Threshold of 1 - kappa", ylab = "MCC", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "MCC")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "MCC"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$F1, xlab = "Threshold of 1 - kappa", ylab = "F1", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "F1")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "F1"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$TPR, xlab = "Threshold of 1 - kappa", ylab = "TPR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "TPR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "TPR"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$FDR, xlab = "Threshold of 1 - kappa", ylab = "FDR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "FDR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "FDR"], lty = 2, col = "blue")
grid()
dev.off()

# Create plot for worse case example when tau tuned.
data_nro <- 5

map <- GHS_MAP_estimation(sim_bdgraph_sf[[data_nro]]$data, verbose = 1, max_iterations = 500,
                          p0 = 50)

kappa <- 1 - map$Kappa

thresholds <- seq(0.25, 0.75, 0.01)
threshold_scores <- data.frame()
for(thr in thresholds) {
  theta_est <- kappa
  theta_est[theta_est < thr] <- 0
  theta_est[theta_est != 0] <- 1
  diag(theta_est) <- 0
  scores <- calculate_scores(conf_matrix(sim_bdgraph_sf[[data_nro]]$theta, theta_est))
  threshold_scores <- rbind(threshold_scores, cbind(thr, scores))
}

pdf("Figures/Supplementary/Thresholding_tau_tuned.pdf", width = 12, height = 3)
par(mfrow = c(1,4))
par(mgp = c(2, 1, 0))
par(mar = c(3.1, 3.1, 0.2, 0.2))
plot(threshold_scores$thr, threshold_scores$MCC, xlab = "Threshold of 1 - kappa", ylab = "MCC", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "MCC")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "MCC"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$F1, xlab = "Threshold of 1 - kappa", ylab = "F1", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "F1")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "F1"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$TPR, xlab = "Threshold of 1 - kappa", ylab = "TPR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "TPR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "TPR"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$FDR, xlab = "Threshold of 1 - kappa", ylab = "FDR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "FDR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "FDR"], lty = 2, col = "blue")
grid()
dev.off()
