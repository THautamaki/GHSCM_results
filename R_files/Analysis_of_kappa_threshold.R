library(GHSGEM)

structures <- c("random", "bdgraph_sf", "huge_sf", "hubs")

# Analysis for the datasets with 100 variables.
p <- 100

# Change "path\\to\\data" part where you have stored datasets.
path <- paste0("path\\to\\data\\p", p, "\\")

path <- "C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\data\\Artikkeli\\p100\\"

load(file = paste0(path, "bdgraph_random_", p, ".Rda"))
load(file = paste0(path, "bdgraph_scale-free_", p, ".Rda"))
load(file = paste0(path, "huge_scale-free_", p, ".Rda"))
load(file = paste0(path, "huge_hubs_", p, ".Rda"))

#######
data_nro <- 2

map <- GHS_MAP_estimate(sim_bdgraph_sf[[data_nro]]$data, verbose = 1, max_iterations = 500)

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

par(mfrow = c(2,2))
par(mgp = c(2, 1, 0))
par(mar = c(3.1, 3.1, 0.2, 0.2))
plot(threshold_scores$thr, threshold_scores$MCC, xlab = "1 - kappa threshold", ylab = "MCC", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "MCC")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "MCC"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$F1, xlab = "1 - kappa threshold", ylab = "F1", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "F1")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "F1"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$TPR, xlab = "1 - kappa threshold", ylab = "TPR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "TPR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "TPR"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$FDR, xlab = "1 - kappa threshold", ylab = "FDR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "FDR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "FDR"], lty = 2, col = "blue")
grid()


#######
data_nro <- 5

map <- GHS_MAP_estimate(sim_bdgraph_sf[[data_nro]]$data, verbose = 1, max_iterations = 500)

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

par(mfrow = c(2,2))
par(mgp = c(2, 1, 0))
par(mar = c(3.1, 3.1, 0.2, 0.2))
plot(threshold_scores$thr, threshold_scores$MCC, xlab = "1 - kappa threshold", ylab = "MCC", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "MCC")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "MCC"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$F1, xlab = "1 - kappa threshold", ylab = "F1", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "F1")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "F1"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$TPR, xlab = "1 - kappa threshold", ylab = "TPR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "TPR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "TPR"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$FDR, xlab = "1 - kappa threshold", ylab = "FDR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "FDR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "FDR"], lty = 2, col = "blue")
grid()


#######
data_nro <- 5

map <- GHS_MAP_estimate(sim_bdgraph_sf[[data_nro]]$data, verbose = 1, max_iterations = 500,
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

par(mfrow = c(2,2))
par(mgp = c(2, 1, 0))
par(mar = c(3.1, 3.1, 0.2, 0.2))
plot(threshold_scores$thr, threshold_scores$MCC, xlab = "1 - kappa threshold", ylab = "MCC", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "MCC")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "MCC"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$F1, xlab = "1 - kappa threshold", ylab = "F1", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "F1")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "F1"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$TPR, xlab = "1 - kappa threshold", ylab = "TPR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "TPR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "TPR"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$FDR, xlab = "1 - kappa threshold", ylab = "FDR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "FDR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "FDR"], lty = 2, col = "blue")
grid()

