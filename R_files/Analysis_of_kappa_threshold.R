library(GHSCM)
library(doParallel)

calculation_over_kappa_thr <- function(p, n = 0, kappa_thresholds = seq(0.25, 0.75, 0.01),
                                       structure = "scale-free", net_size = 0,
                                       n_repeats = 10, n_cores = 0, seed = 22102025) {
  if (n == 0) n <- p * 0.8
  if (n_cores == 0) n_cores <- min(c(n_repeats, detectCores() - 1))
  sims <- list()
  set.seed(seed)
  for (r in 1:n_repeats) {
    if (structure == "scale-free" | structure == "hub" | structure == "cluster") {
      sim <- huge::huge.generator(n = n, d = p, graph = structure, verbose = FALSE)
    }
    else if (structure == "random") {
      if (net_size == 0) net_size <- floor(p/2)
      sim <- BDgraph::bdgraph.sim(p = p, n = n, graph = "random", size = net_size)
      sim$omega <- sim$K
      sim$theta <- as.matrix(sim$G)
      sim$sparsity <- (sum(sim$theta) / 2) / (p * (p - 1) / 2)
    }
    sims[[r]] <- sim
  }
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  all_scores <- foreach (r = 1:n_repeats, .packages = c("GHSCM", "huge", "BDgraph"), .combine = rbind) %dopar% {
    sim <- sims[[r]]
    map <- GHSCM::GHS_MAP_estimation(sim$data)
    true_sparsity <- 1 - sim$sparsity
    kappa <- 1 - map$Kappa
    one_round_scores <- data.frame()
    for (thr in kappa_thresholds) {
      theta_est <- kappa
      theta_est[theta_est < thr] <- 0
      theta_est[theta_est != 0] <- 1
      diag(theta_est) <- 0
      sparsity <- 1 - (sum(theta_est) / 2) / (p * (p - 1)/2)
      scores <- calculate_scores(conf_matrix(sim$theta, theta_est))
      one_round_scores <- rbind(one_round_scores, cbind(r, thr, scores, sparsity, true_sparsity))
    }
    one_round_scores
  }
  stopCluster(cl)
  stopImplicitCluster()
  return(all_scores)
}

plot_sparsity_mcc_f1 <- function(p, mean_scores, true_sparsity, x_label, ylim_sparsity, legend,
                                 legend_position = "bottomleft", legend_in_plot = 1,
                                 plot_rows = 1, filename = NULL, pdf_size = c(10.5, 3.5)) {
  if (!is.null(filename)) pdf(filename, width = pdf_size[1], height = pdf_size[2])
  ltys <- rep(c(2,4,5,6,1), length.out = length(mean_scores))
  cols <- palette.colors(length(mean_scores))
  if (plot_rows == 1) par(mfrow = c(1,3))
  par(mgp = c(2, 1, 0))
  par(mar = c(3.1, 3.1, 0.5, 0.2))
  plot(mean_scores[[1]]$thr, mean_scores[[1]]$sparsity, type = "l", xlab = x_label,
       ylab = "sparsity", lty = ltys[1], col = cols[1], ylim = ylim_sparsity)
  for (i in 2:length(mean_scores)) {
    lines(mean_scores[[i]]$thr, mean_scores[[i]]$sparsity, lty = ltys[i], col = cols[i])
  }
  grid()
  abline(h = true_sparsity, lty = 2, col = "blue")
  abline(v = 0.5, lty = 2, col = "blue")
  if (legend_in_plot == 1) {
    legend(legend_position, legend, lty = c(ltys, 2), col = c(cols, "blue"), bg = "white")
  }
  
  plot(mean_scores[[1]]$thr, mean_scores[[1]]$MCC, type = "l", xlab = x_label, ylab = "MCC",
       ylim = c(0,1), lty = ltys[1], col = cols[1])
  for (i in 2:length(mean_scores)) {
    lines(mean_scores[[i]]$thr, mean_scores[[i]]$MCC, lty = ltys[i], col = cols[i])
  }
  grid()
  abline(v = 0.5, lty = 2, col = "blue")
  if (legend_in_plot == 2) {
    legend(legend_position, legend, lty = c(ltys, 2), col = c(cols, "blue"), bg = "white")
  }
  
  plot(mean_scores[[1]]$thr, mean_scores[[1]]$F1, type = "l", xlab = x_label,
       ylab = expression(paste(F[1], "-score")), ylim = c(0,1), lty = ltys[1], col = cols[1])
  for (i in 2:length(mean_scores)) {
    lines(mean_scores[[i]]$thr, mean_scores[[i]]$F1, lty = ltys[i], col = cols[i])
  }
  grid()
  abline(v = 0.5, lty = 2, col = "blue")
  if (legend_in_plot == 3) {
    legend(legend_position, legend, lty = c(ltys, 2), col = c(cols, "blue"), bg = "white")
  }
  if (!is.null(filename)) dev.off()
}

plot_mcc_tpr_fdr <- function(p, mean_scores, true_sparsity, x_label, legend,
                             legend_position = "bottomleft", legend_in_plot = 1,
                             plot_rows = 1, filename = NULL, save_eps = FALSE, save_pdf = FALSE,
                             image_size = c(10.5, 3.5)) {
  if (save_eps) {
    setEPS()
    postscript(filename, width = image_size[1], height = image_size[2])
  }
  else if (save_pdf) {
    pdf(filename, width = image_size[1], height = image_size[2])
  }
  ltys <- rep(c(2,4,5,6,1), length.out = length(mean_scores))
  cols <- palette.colors(length(mean_scores))
  if (plot_rows == 1) par(mfrow = c(1,3))
  par(mgp = c(2, 1, 0))
  par(mar = c(3.1, 3.3, 0.5, 0.2))
  plot(mean_scores[[1]]$thr, mean_scores[[1]]$MCC, type = "l", xlab = x_label,
       ylab = "MCC", lty = ltys[1], col = cols[1], ylim = c(0,1))
  for (i in 2:length(mean_scores)) {
    lines(mean_scores[[i]]$thr, mean_scores[[i]]$MCC, lty = ltys[i], col = cols[i])
  }
  grid()
  abline(v = p - 1, lty = 2, col = "blue")
  if (legend_in_plot == 1) {
    legend(legend_position, legend, lty = c(ltys, 2), col = c(cols, "blue"), bg = "white")
  }
  
  plot(mean_scores[[1]]$thr, mean_scores[[1]]$TPR, type = "l", xlab = x_label, ylab = "TPR",
       ylim = c(0,1), lty = ltys[1], col = cols[1])
  for (i in 2:length(mean_scores)) {
    lines(mean_scores[[i]]$thr, mean_scores[[i]]$TPR, lty = ltys[i], col = cols[i])
  }
  grid()
  abline(v = p - 1, lty = 2, col = "blue")
  if (legend_in_plot == 2) {
    legend(legend_position, legend, lty = c(ltys, 2), col = c(cols, "blue"), bg = "white")
  }
  
  plot(mean_scores[[1]]$thr, mean_scores[[1]]$FDR, type = "l", xlab = x_label, ylab = "FDR",
       ylim = c(0,1), lty = ltys[1], col = cols[1])
  for (i in 2:length(mean_scores)) {
    lines(mean_scores[[i]]$thr, mean_scores[[i]]$FDR, lty = ltys[i], col = cols[i])
  }
  grid()
  abline(v = p - 1, lty = 2, col = "blue")
  if (legend_in_plot == 2) {
    legend(legend_position, legend, lty = c(ltys, 2), col = c(cols, "blue"), bg = "white")
  }
  if (!is.null(filename)) dev.off()
}

# Set data dimensions.
n <- 120
p <- 100

# Set path (no needed to change).
path <- paste0("Data/n", n, "_p", p, "/")

# Load datasets.
sim_bdgraph_sf <- readRDS(file = paste0(path, "bdgraph_scale-free_n", n, "_p", p, ".Rds"))

# Create plot for normal case example.
data_nro <- 2

map <- GHS_MAP_estimation(sim_bdgraph_sf[[data_nro]]$data, verbose = 0)

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

x_label <- expression(paste("Threshold of 1 - ", kappa))

pdf("Figures/Supplementary/Thresholding_normal_case_example.pdf", width = 12, height = 3)
par(mfrow = c(1,4))
par(mgp = c(2, 1, 0))
par(mar = c(3.1, 3.1, 0.2, 0.2))
plot(threshold_scores$thr, threshold_scores$MCC, xlab = x_label, ylab = "MCC", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "MCC")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "MCC"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$F1, xlab = x_label, ylab = "F1", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "F1")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "F1"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$TPR, xlab = x_label, ylab = "TPR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "TPR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "TPR"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$FDR, xlab = x_label, ylab = "FDR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "FDR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "FDR"], lty = 2, col = "blue")
grid()
dev.off()

# Create plot for worse case example.
data_nro <- 5

map <- GHS_MAP_estimation(sim_bdgraph_sf[[data_nro]]$data, verbose = 0)

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
plot(threshold_scores$thr, threshold_scores$MCC, xlab = x_label, ylab = "MCC", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "MCC")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "MCC"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$F1, xlab = x_label, ylab = "F1", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "F1")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "F1"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$TPR, xlab = x_label, ylab = "TPR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "TPR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "TPR"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$FDR, xlab = x_label, ylab = "FDR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "FDR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "FDR"], lty = 2, col = "blue")
grid()
dev.off()

# Create plot for worse case example when tau tuned.
data_nro <- 5

map <- GHS_MAP_estimation(sim_bdgraph_sf[[data_nro]]$data, verbose = 0, p0 = 50)

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
plot(threshold_scores$thr, threshold_scores$MCC, xlab = x_label, ylab = "MCC", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "MCC")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "MCC"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$F1, xlab = x_label, ylab = "F1", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "F1")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "F1"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$TPR, xlab = x_label, ylab = "TPR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "TPR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "TPR"], lty = 2, col = "blue")
grid()
plot(threshold_scores$thr, threshold_scores$FDR, xlab = x_label, ylab = "FDR", type = "l")
points(threshold_scores[threshold_scores$thr == 0.5, c("thr", "FDR")], cex = 2, col = "blue")
abline(v = 0.5, lty = 2, col = "blue")
abline(h = threshold_scores[threshold_scores$thr == 0.5, "FDR"], lty = 2, col = "blue")
grid()
dev.off()


#####
### Analysis with p = 100 and scale-free structure.
# Define needed parameters.
p <- 100
ratios <- c(0.6, 0.8, 1, 1.25, 1.5, 2, 4)
n_repeats <- 20
structure = "scale-free"
seed <- 20251211

# Run analysis if results file does not exist.
if (!file.exists("Results_files/Kappa_analyses/kappa_scores_p100_sf.rds")) {
  print(s1 <- Sys.time())
  scores_p100_sf <- list()
  for (ratio in ratios) {
    n <- floor(p * ratio)
    cat("n =", n, "\n")
    scores_p100_sf[[paste0("n", n)]] <- calculation_over_kappa_thr(p = p, n = n,
                                                                   n_repeats = n_repeats,
                                                                   structure = structure,
                                                                   seed = seed)
  }
  print(Sys.time() - s1)
  saveRDS(kappa_scores_p100_sf, file = "Results_files/Kappa_analyses/kappa_scores_p100_sf.rds")
}

# Load results if analysis is already done.
if (file.exists("Results_files/Kappa_analyses/kappa_scores_p100_sf.rds")) {
  scores_p100_sf <- readRDS("Results_files/Kappa_analyses/kappa_scores_p100_sf.rds")
}

# Calculate mean scores.
mean_scores_p100_sf <- list()
for (ratio in ratios) {
  n <- floor(p * ratio)
  mean_scores_p100_sf[[paste0("n", n)]] <- aggregate(kappa_scores_p100_sf[[paste0("n", n)]],
                                                     by = list(kappa_scores_p100_sf[[paste0("n", n)]]$thr),
                                                     FUN = mean)
}

# Calculate true sparsity and combine estimated sparsities for y-axis limits.
true_sparsities <- est_sparsities_sf_p100 <- data.frame()
for (i in 1:length(mean_scores_p100_sf)) {
  true_sparsities <- rbind(true_sparsities, mean_scores_p100_sf[[i]]$true_sparsity)
  est_sparsities_sf_p100 <- rbind(est_sparsities_sf_p100, mean_scores_p100_sf[[i]]$sparsity)
}
if (all(true_sparsities == true_sparsities[1,1])) {
  true_sparsity_sf_p100 <- true_sparsities[1,1]
  cat("All sparsities are equal. True sparsity is:", true_sparsity_sf_p100)
} else {
  true_sparsity_sf_p100 <- mean(as.matrix(true_sparsities))
  sd_sparsity <- sd(as.matrix(true_sparsities))
  cat("Sparsities are not equal. Mean (sd) of true sparsity is: ", true_sparsity_sf_p100,
      " (", round(sd_sparsity, 6), ")", sep = "")
}

# Define legend for the plot and label for x-axis.
legend_p100 <- c(paste("n =", floor(p * ratios)), "true sparsity / default thr.")
x_label <- expression(paste("Threshold of 1 - ", kappa))

# Create plot.
plot_sparsity_mcc_f1(p, mean_scores_p100_sf, true_sparsity = true_sparsity_sf_p100, x_label = x_label,
                     legend = legend_p100, ylim_sparsity = c(min(est_sparsities_sf_p100), max(est_sparsities_sf_p100)))

### Analysis with p = 100 and random structure with network size p/2.
# Define needed parameters.
p <- 100
ratios <- c(0.6, 0.8, 1, 1.25, 1.5, 2, 4)
n_repeats <- 20
structure = "random"
seed <- 20251212

# Run analysis if results file does not exist.
if (!file.exists("Results_files/Kappa_analyses/kappa_scores_p100_random.rds")) {
  print(s2 <- Sys.time())
  scores_p100_random <- list()
  for (ratio in ratios) {
    n <- floor(p * ratio)
    cat("n =", n, "\n")
    scores_p100_random[[paste0("n", n)]] <- calculation_over_kappa_thr(p = p, n = n,
                                                                       structure = structure,
                                                                       seed = seed)
  }
  print(Sys.time() - s2)
  saveRDS(scores_p100_random, file = "Results_files/Kappa_analyses/kappa_scores_p100_random.rds")
}

# Load results if analysis is already done.
if (file.exists("Results_files/Kappa_analyses/kappa_scores_p100_random.rds")) {
  scores_p100_random <- readRDS("Results_files/Kappa_analyses/kappa_scores_p100_random.rds")
}

# Calculate mean scores.
mean_scores_p100_random <- list()
for (ratio in ratios) {
  n <- floor(p * ratio)
  mean_scores_p100_random[[paste0("n", n)]] <- aggregate(scores_p100_random[[paste0("n", n)]],
                                                         by = list(scores_p100_random[[paste0("n", n)]]$thr),
                                                         FUN = mean)
}

# Calculate true sparsity and combine estimated sparsities for y-axis limits.
true_sparsities <- est_sparsities_random_p100 <- data.frame()
for (i in 1:length(mean_scores_p100_random)) {
  true_sparsities <- rbind(true_sparsities, mean_scores_p100_random[[i]]$true_sparsity)
  est_sparsities_random_p100 <- rbind(est_sparsities_random_p100, mean_scores_p100_random[[i]]$sparsity)
}
if (all(true_sparsities == true_sparsities[1,1])) {
  true_sparsity_random_p100 <- true_sparsities[1,1]
  cat("All sparsities are equal. True sparsity is:", true_sparsity_random_p100)
} else {
  true_sparsity_random_p100 <- mean(as.matrix(true_sparsities))
  sd_sparsity <- sd(as.matrix(true_sparsities))
  cat("Sparsities are not equal. Mean (sd) of true sparsity is: ", true_sparsity_random_p100,
      " (", round(sd_sparsity, 6), ")", sep = "")
}

# Define legend for the plot and label for x-axis.
legend_p100 <- c(paste("n =", floor(p * ratios)), "true sparsity / default thr.")
x_label <- expression(paste("Threshold of 1 - ", kappa))

# Create plot.
plot_sparsity_mcc_f1(p, mean_scores_p100_random, true_sparsity = true_sparsity_random_p100, x_label = x_label,
                     legend = legend_p100, ylim_sparsity = c(min(est_sparsities_random_p100), max(est_sparsities_random_p100)))


### Analysis with p = 100 and random structure with network size p * 1.5.
# Define needed parameters.
p <- 100
ratios <- c(0.6, 0.8, 1, 1.25, 1.5, 2, 4)
net_size <- floor(p * 1.5)
structure = "random"
seed <- 20251213

# Run analysis if results file does not exist.
if (!file.exists("Results_files/Kappa_analyses/kappa_scores_p100_random_denser.rds")) {
  print(s3 <- Sys.time())
  scores_p100_random_denser <- list()
  for (ratio in ratios) {
    n <- floor(p * ratio)
    cat("n =", n, "\n")
    scores_p100_random_denser[[paste0("n", n)]] <- calculation_over_kappa_thr(p = p, n = n,
                                                                              n_repeats = n_repeats,
                                                                              structure = structure,
                                                                              net_size = net_size,
                                                                              seed = seed)
  }
  print(Sys.time() - s3)
  saveRDS(scores_p100_random_denser, file = "Results_files/Kappa_analyses/kappa_scores_p100_random_denser.rds")
}

# Load results if analysis is already done.
if (file.exists("Results_files/Kappa_analyses/kappa_scores_p100_random_denser.rds")) {
  scores_p100_random_denser <- readRDS("Results_files/Kappa_analyses/kappa_scores_p100_random_denser.rds")
}

# Calculate mean scores.
mean_scores_p100_random_denser <- list()
for (ratio in ratios) {
  n <- floor(p * ratio)
  mean_scores_p100_random_denser[[paste0("n", n)]] <- aggregate(scores_p100_random_denser[[paste0("n", n)]],
                                                                by = list(scores_p100_random_denser[[paste0("n", n)]]$thr),
                                                                FUN = mean)
}

# Calculate true sparsity and combine estimated sparsities for y-axis limits.
true_sparsities <- est_sparsities_rnd_dens_p100 <- data.frame()
for (i in 1:length(mean_scores_p100_random_denser)) {
  true_sparsities <- rbind(true_sparsities, mean_scores_p100_random_denser[[i]]$true_sparsity)
  est_sparsities_rnd_dens_p100 <- rbind(est_sparsities_rnd_dens_p100, mean_scores_p100_random_denser[[i]]$sparsity)
}
if (all(true_sparsities == true_sparsities[1,1])) {
  true_sparsity_rnd_dens_p100 <- true_sparsities[1,1]
  cat("All sparsities are equal. True sparsity is:", true_sparsity_rnd_dens_p100)
} else {
  true_sparsity_rnd_dens_p100 <- mean(as.matrix(true_sparsities))
  sd_sparsity <- sd(as.matrix(true_sparsities))
  cat("Sparsities are not equal. Mean (sd) of true sparsity is: ", true_sparsity_rnd_dens_p100,
      " (", round(sd_sparsity, 6), ")", sep = "")
}

# Define legend for the plot and label for x-axis.
legend_p100 <- c(paste("n =", floor(p * ratios)), "true sparsity / default thr.")
x_label <- expression(paste("Threshold of 1 - ", kappa))

# Create plot.
plot_sparsity_mcc_f1(p, mean_scores_p100_random_denser, true_sparsity = true_sparsity_rnd_dens_p100, x_label = x_label,
                     legend = legend_p100, ylim_sparsity = c(min(est_sparsities_rnd_dens_p100), max(est_sparsities_rnd_dens_p100)))


### Analysis with p = 100 and hub structure.
# Define needed parameters.
p <- 100
ratios <- c(0.6, 0.8, 1, 1.25, 1.5, 2, 4)
n_repeats <- 20
structure = "hub"
seed <- 20251214

# Run analysis if results file does not exist.
if (!file.exists("Results_files/Kappa_analyses/kappa_scores_p100_hub.rds")) {
  print(s6 <- Sys.time())
  scores_p100_hub <- list()
  for (ratio in ratios) {
    n <- floor(p * ratio)
    cat("n =", n, "\n")
    scores_p100_hub[[paste0("n", n)]] <- calculation_over_kappa_thr(p = p, n = n,
                                                                    n_repeats = n_repeats,
                                                                    structure = structure,
                                                                    seed = seed)
  }
  print(Sys.time() - s6)
  saveRDS(scores_p100_hub, file = "Results_files/Kappa_analyses/kappa_scores_p100_hub.rds")
}

# Load results if analysis is already done.
if (file.exists("Results_files/Kappa_analyses/kappa_scores_p100_hub.rds")) {
  scores_p100_hub <- readRDS("Results_files/Kappa_analyses/kappa_scores_p100_hub.rds")
}

# Calculate mean scores.
mean_scores_p100_hub <- list()
for (ratio in ratios) {
  n <- floor(p * ratio)
  mean_scores_p100_hub[[paste0("n", n)]] <- aggregate(scores_p100_hub[[paste0("n", n)]],
                                                      by = list(scores_p100_hub[[paste0("n", n)]]$thr),
                                                      FUN = mean)
}

# Calculate true sparsity and combine estimated sparsities for y-axis limits.
true_sparsities <- est_sparsities_hub_p100 <- data.frame()
for (i in 1:length(mean_scores_p100_hub)) {
  true_sparsities <- rbind(true_sparsities, mean_scores_p100_hub[[i]]$true_sparsity)
  est_sparsities_hub_p100 <- rbind(est_sparsities_hub_p100, mean_scores_p100_hub[[i]]$sparsity)
}
if (all(true_sparsities == true_sparsities[1,1])) {
  true_sparsity_hub_p100 <- true_sparsities[1,1]
  cat("All sparsities are equal. True sparsity is:", true_sparsity_hub_p100)
} else {
  true_sparsity_hub_p100 <- mean(as.matrix(true_sparsities))
  sd_sparsity <- sd(as.matrix(true_sparsities))
  cat("Sparsities are not equal. Mean (sd) of true sparsity is: ", true_sparsity_hub_p100,
      " (", round(sd_sparsity, 6), ")", sep = "")
}

# Define legend for the plot and label for x-axis.
legend_p100 <- c(paste("n =", floor(p * ratios)), "true sparsity / default thr.")
x_label <- expression(paste("Threshold of 1 - ", kappa))

# Create plot.
plot_sparsity_mcc_f1(p, mean_scores_p100_hub, true_sparsity = true_sparsity_hub_p100, x_label = x_label,
                     legend = legend_p100, ylim_sparsity = c(min(est_sparsities_hub_p100), max(est_sparsities_hub_p100)))

## Combine plots with p = 100 and scale-free and random structures in one figure.
pdf("Figures/Supplementary/Kappa_analysis_p100.pdf", width = 10.5, height = 10.5)
par(mfrow = c(3,3))
plot_sparsity_mcc_f1(100, mean_scores_p100_sf, true_sparsity = true_sparsity_sf_p100,
                     x_label = x_label, legend = legend_p100, plot_rows = 3, legend_position = "bottomright",
                     ylim_sparsity = c(min(est_sparsities_sf_p100), max(est_sparsities_sf_p100)))
plot_sparsity_mcc_f1(100, mean_scores_p100_random, true_sparsity = true_sparsity_random_p100,
                     x_label = x_label, legend = legend_p100, plot_rows = 3, legend_position = "bottomright",
                     ylim_sparsity = c(min(est_sparsities_random_p100), max(est_sparsities_random_p100)))
plot_sparsity_mcc_f1(100, mean_scores_p100_random_denser, true_sparsity = true_sparsity_rnd_dens_p100,
                     x_label = x_label, legend = legend_p100, plot_rows = 3, legend_position = "bottomright",
                     ylim_sparsity = c(min(est_sparsities_rnd_dens_p100), max(est_sparsities_rnd_dens_p100)))

dev.off()


### Analysis with p = 200 and hub structure.
# Define needed parameters.
p <- 200
ratios <- c(0.6, 0.8, 1, 1.25, 1.5, 2, 4)
structure = "hub"
seed <- 20251215

# Run analysis if results file does not exist.
if (!file.exists("Results_files/Kappa_analyses/kappa_scores_p200_hub.rds")) {
  print(s4 <- Sys.time())
  scores_p200_hub <- list()
  for (ratio in ratios) {
    n <- floor(p * ratio)
    cat("n =", n, "\n")
    scores_p200_hub[[paste0("n", n)]] <- calculation_over_kappa_thr(p = p, n = n,
                                                                    n_repeats = n_repeats,
                                                                    structure = structure,
                                                                    seed = seed)
  }
  print(Sys.time() - s4)
  saveRDS(scores_p200_hub, file = "Results_files/Kappa_analyses/kappa_scores_p200_hub.rds")
}

# Load results if analysis is already done.
if (file.exists("Results_files/Kappa_analyses/kappa_scores_p200_hub.rds")) {
  scores_p200_hub <- readRDS("Results_files/Kappa_analyses/kappa_scores_p200_hub.rds")
}

# Calculate mean scores.
mean_scores_p200_hub <- list()
for (ratio in ratios) {
  n <- floor(p * ratio)
  mean_scores_p200_hub[[paste0("n", n)]] <- aggregate(scores_p200_hub[[paste0("n", n)]],
                                                      by = list(scores_p200_hub[[paste0("n", n)]]$thr),
                                                      FUN = mean)
}

# Calculate true sparsity and combine estimated sparsities for y-axis limits.
true_sparsities <- est_sparsities_hub_p200 <- data.frame()
for (i in 1:length(mean_scores_p200_hub)) {
  true_sparsities <- rbind(true_sparsities, mean_scores_p200_hub[[i]]$true_sparsity)
  est_sparsities_hub_p200 <- rbind(est_sparsities_hub_p200, mean_scores_p200_hub[[i]]$sparsity)
}
if (all(true_sparsities == true_sparsities[1,1])) {
  true_sparsity_hub_p200 <- true_sparsities[1,1]
  cat("All sparsities are equal. True sparsity is:", true_sparsity_hub_p200)
} else {
  true_sparsity_hub_p200 <- mean(as.matrix(true_sparsities))
  sd_sparsity <- sd(as.matrix(true_sparsities))
  cat("Sparsities are not equal. Mean (sd) of true sparsity is: ", true_sparsity_hub_p200,
      " (", round(sd_sparsity, 6), ")", sep = "")
}

# Define legend for the plot and label for x-axis.
legend_p200 <- c(paste("n =", floor(p * ratios)), "true sparsity / default thr.") 
x_label <- expression(paste("Threshold of 1 - ", kappa))

# Create and save plot.
plot_sparsity_mcc_f1(p, mean_scores_p200_hub, true_sparsity = true_sparsity_hub_p200, x_label = x_label,
                     legend = legend_p200, legend_position = "topleft",
                     ylim_sparsity = c(min(est_sparsities_hub_p200), max(est_sparsities_hub_p200)))


### Analysis with p = 250 and cluster structure.
# Define needed parameters.
p <- 250
ratios <- c(0.6, 0.8, 1, 1.25, 1.5, 2, 4)
n_repeats <- 10
structure = "cluster"
seed <- 20251216

# Run analysis if results file does not exist.
if (!file.exists("Results_files/Kappa_analyses/kappa_scores_p250_cluster.rds")) {
  print(s5 <- Sys.time())
  scores_p250_cluster <- list()
  for (ratio in ratios) {
    n <- floor(p * ratio)
    cat("n =", n, "\n")
    scores_p250_cluster[[paste0("n", n)]] <- calculation_over_kappa_thr(p = p, n = n,
                                                                        n_repeats = n_repeats,
                                                                        structure = structure,
                                                                        seed = seed)
  }
  print(Sys.time() - s5)
  saveRDS(scores_p250_cluster, file = "Results_files/Kappa_analyses/kappa_scores_p250_cluster.rds")
}

# Load results if analysis is already done.
if (file.exists("Results_files/Kappa_analyses/kappa_scores_p250_cluster.rds")) {
  scores_p250_cluster <- readRDS("Results_files/Kappa_analyses/kappa_scores_p250_cluster.rds")
}

# Calculate mean scores.
mean_scores_p250_cluster <- list()
for (ratio in ratios) {
  n <- floor(p * ratio)
  mean_scores_p250_cluster[[paste0("n", n)]] <- aggregate(scores_p250_cluster[[paste0("n", n)]],
                                                          by = list(scores_p250_cluster[[paste0("n", n)]]$thr),
                                                          FUN = mean)
}

# Calculate true sparsity and combine estimated sparsities for y-axis limits.
true_sparsities <- est_sparsities_cluster_p250 <- data.frame()
for (i in 1:length(mean_scores_p250_cluster)) {
  true_sparsities <- rbind(true_sparsities, mean_scores_p250_cluster[[i]]$true_sparsity)
  est_sparsities_cluster_p250 <- rbind(est_sparsities_cluster_p250, mean_scores_p250_cluster[[i]]$sparsity)
}
if (all(true_sparsities == true_sparsities[1,1])) {
  true_sparsity_cluster_p250 <- true_sparsities[1,1]
  cat("All sparsities are equal. True sparsity is:", true_sparsity_cluster_p250)
} else {
  true_sparsity_cluster_p250 <- mean(as.matrix(true_sparsities))
  sd_sparsity <- sd(as.matrix(true_sparsities))
  cat("Sparsities are not equal. Mean (sd) of true sparsity is: ", true_sparsity_cluster_p250, " (", round(sd_sparsity, 6), ")", sep = "")
}

# Define legend for the plot and label for x-axis.
legend_p250 <- c(paste("n =", floor(p * ratios)), "true sparsity / default thr.") 
x_label <- expression(paste("Threshold of 1 - ", kappa))

# Create and save plot.
plot_sparsity_mcc_f1(p, mean_scores_p250_cluster, true_sparsity = true_sparsity_cluster_p250, x_label = x_label,
                     legend = legend_p200, legend_position = "bottomleft", legend_in_plot = 2,
                     ylim_sparsity = c(min(est_sparsities_cluster_p250), max(est_sparsities_cluster_p250)))


# Create plot.
## Combine plots with p = 100 and scale-free and random structures in one figure.
pdf("Figures/Supplementary/Kappa_analysis_p100-200_hub_p250_cluster.pdf", width = 10.5, height = 10.5)
par(mfrow = c(3,3))
plot_sparsity_mcc_f1(100, mean_scores_p100_hub, true_sparsity = true_sparsity_hub_p100, x_label = x_label,
                     legend = legend_p100, ylim_sparsity = c(min(est_sparsities_hub_p100), max(est_sparsities_hub_p100)),
                     plot_rows = 3, legend_position = "bottomright")

plot_sparsity_mcc_f1(200, mean_scores_p200_hub, true_sparsity = true_sparsity_hub_p200, x_label = x_label,
                     legend = legend_p200, legend_position = "bottomright",
                     ylim_sparsity = c(min(est_sparsities_hub_p200), max(est_sparsities_hub_p200)),
                     plot_rows = 3)

plot_sparsity_mcc_f1(250, mean_scores_p250_cluster, true_sparsity = true_sparsity_cluster_p250, x_label = x_label,
                     legend = legend_p250, legend_position = "bottomleft", legend_in_plot = 2,
                     ylim_sparsity = c(min(est_sparsities_cluster_p250), max(est_sparsities_cluster_p250)),
                     plot_rows = 3)
dev.off()
