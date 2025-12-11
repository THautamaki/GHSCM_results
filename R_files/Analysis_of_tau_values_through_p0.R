library(doParallel)
library(latex2exp)

calculation_over_p0 <- function(p, n = 0, p0s = NULL, n_p0s = 50, step_p0 = 5,
                                structure = "scale-free", net_size = 0,
                                n_repeats = 10, n_cores = 0, seed = 22102025) {
  if (n == 0) n <- p * 0.8
  if (is.null(p0s)) p0s <- seq(ifelse(p - 1 > n_p0s, p - 1 - n_p0s, 1), p - 1 + n_p0s, step_p0)
  if (n_cores == 0) n_cores <- min(c(length(p0s), detectCores() - 1))
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
  all_scores <- foreach (p0 = p0s, .packages = c("GHSCM", "huge", "BDgraph"), .combine = rbind) %dopar% {
    one_round_scores <- data.frame()
    for (r in 1:n_repeats) {
      sim <- sims[[r]]
      data <- sim$data
      map <- GHSCM::GHS_MAP_estimation(data, p0 = p0)
      sparsity <- 1 - (sum(map$Theta_est) / 2) / (p * (p - 1)/2)
      cm <- conf_matrix(sim$theta, map$Theta_est)
      scores <- calculate_scores(cm)
      sl_omega <- stein_loss(sim$omega, map$Omega_est)
      sl_sigma <- stein_loss(sim$sigma, map$Sigma_est)
      f_norm_omega <- norm(sim$omega - map$Omega_est, type = "f")
      f_norm_sigma <- norm(sim$sigma - map$Sigma_est, type = "f")
      f_norm_rel <- f_norm_omega / norm(sim$omega, type = "f")
      iters <- map$iters
      true_sparsity <- 1 - sim$sparsity
      one_round_scores <- rbind(one_round_scores, cbind(r, p0, scores, sparsity, sl_omega, f_norm_omega,
                                                        f_norm_rel, sl_sigma, f_norm_sigma, iters,
                                                        true_sparsity))
    }
    one_round_scores
  }
  stopCluster(cl)
  stopImplicitCluster()
  return(all_scores)
}

calculate_tau <- function(n, p, p0) {
  c <- p/5
  if (p < 250) {
    c <- 50
  }
  (p0/(p * (p - 1)/2)) * (c * sqrt(p0)/(n * sqrt(n)))
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
  plot(mean_scores[[1]]$p0, mean_scores[[1]]$sparsity, type = "l", xlab = x_label,
       ylab = "sparsity", lty = ltys[1], col = cols[1], ylim = ylim_sparsity)
  for (i in 2:length(mean_scores)) {
    lines(mean_scores[[i]]$p0, mean_scores[[i]]$sparsity, lty = ltys[i], col = cols[i])
  }
  grid()
  abline(h = true_sparsity, lty = 2, col = "blue")
  abline(v = p - 1, lty = 2, col = "blue")
  if (legend_in_plot == 1) {
    legend(legend_position, legend, lty = c(ltys, 2), col = c(cols, "blue"), bg = "white")
  }
  
  plot(mean_scores[[1]]$p0, mean_scores[[1]]$MCC, type = "l", xlab = x_label, ylab = "MCC",
       ylim = c(0,1), lty = ltys[1], col = cols[1])
  for (i in 2:length(mean_scores)) {
    lines(mean_scores[[i]]$p0, mean_scores[[i]]$MCC, lty = ltys[i], col = cols[i])
  }
  grid()
  abline(v = p - 1, lty = 2, col = "blue")
  if (legend_in_plot == 2) {
    legend(legend_position, legend, lty = c(ltys, 2), col = c(cols, "blue"), bg = "white")
  }
  
  plot(mean_scores[[1]]$p0, mean_scores[[1]]$F1, type = "l", xlab = x_label,
       ylab = expression(paste(F[1], "-score")), ylim = c(0,1), lty = ltys[1], col = cols[1])
  for (i in 2:length(mean_scores)) {
    lines(mean_scores[[i]]$p0, mean_scores[[i]]$F1, lty = ltys[i], col = cols[i])
  }
  grid()
  abline(v = p - 1, lty = 2, col = "blue")
  if (legend_in_plot == 3) {
    legend(legend_position, legend, lty = c(ltys, 2), col = c(cols, "blue"), bg = "white")
  }
  if (!is.null(filename)) dev.off()
}

plot_sparsity_mcc_sl <- function(p, mean_scores, true_sparsity, x_label, ylim_sparsity,
                                 ylim_steinloss, legend,
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
  plot(mean_scores[[1]]$p0, mean_scores[[1]]$sparsity, type = "l", xlab = x_label,
       ylab = "sparsity", lty = ltys[1], col = cols[1], ylim = ylim_sparsity)
  for (i in 2:length(mean_scores)) {
    lines(mean_scores[[i]]$p0, mean_scores[[i]]$sparsity, lty = ltys[i], col = cols[i])
  }
  grid()
  abline(h = true_sparsity, lty = 2, col = "blue")
  abline(v = p - 1, lty = 2, col = "blue")
  if (legend_in_plot == 1) {
    legend(legend_position, legend, lty = c(ltys, 2), col = c(cols, "blue"), bg = "white")
  }
  
  plot(mean_scores[[1]]$p0, mean_scores[[1]]$MCC, type = "l", xlab = x_label, ylab = "MCC",
       ylim = c(0,1), lty = ltys[1], col = cols[1])
  for (i in 2:length(mean_scores)) {
    lines(mean_scores[[i]]$p0, mean_scores[[i]]$MCC, lty = ltys[i], col = cols[i])
  }
  grid()
  abline(v = p - 1, lty = 2, col = "blue")
  if (legend_in_plot == 2) {
    legend(legend_position, legend, lty = c(ltys, 2), col = c(cols, "blue"), bg = "white")
  }
  
  plot(mean_scores[[1]]$p0, mean_scores[[1]]$sl_omega, type = "l", xlab = x_label,
       ylab = TeX("Stein's loss ($\\widehat{\\Omega}$)"), ylim = ylim_steinloss, lty = ltys[1], col = cols[1])
  for (i in 2:length(mean_scores)) {
    lines(mean_scores[[i]]$p0, mean_scores[[i]]$sl_omega, lty = ltys[i], col = cols[i])
  }
  grid()
  abline(v = p - 1, lty = 2, col = "blue")
  if (legend_in_plot == 3) {
    legend(legend_position, legend, lty = c(ltys, 2), col = c(cols, "blue"), bg = "white")
  }
  if (!is.null(filename)) dev.off()
}

### Analysis with p = 100 and scale-free structure.
# Define needed parameters.
p <- 100
ratios <- c(0.6, 0.8, 1, 1.25, 1.5, 2, 4)
n_repeats <- 20
n_p0s <- 75
step_p0 <- 5
structure = "scale-free"
seed <- 22102025

# Run analysis if results file does not exist.
if (!file.exists("Results_files/Tau_analyses/scores_p100_sf.rds")) {
  print(s1 <- Sys.time())
  scores_p100_sf <- list()
  for (ratio in ratios) {
    n <- floor(p * ratio)
    cat("n =", n, "\n")
    scores_p100_sf[[paste0("n", n)]] <- calculation_over_p0(p = p, n = n,
                                                            n_repeats = n_repeats,
                                                            n_p0s = n_p0s,
                                                            step_p0 = step_p0,
                                                            structure = structure,
                                                            seed = seed)
  }
  print(Sys.time() - s1)
  saveRDS(scores_p100_sf, file = "Results_files/Tau_analyses/scores_p100_sf.rds")
}

# Load results if analysis is already done.
if (file.exists("Results_files/Tau_analyses/scores_p100_sf.rds")) {
  scores_p100_sf <- readRDS("Results_files/Tau_analyses/scores_p100_sf.rds")
}

# Calculate mean scores.
mean_scores_p100_sf <- list()
for (ratio in ratios) {
  n <- floor(p * ratio)
  mean_scores_p100_sf[[paste0("n", n)]] <- aggregate(scores_p100_sf[[paste0("n", n)]],
                                                     by = list(scores_p100_sf[[paste0("n", n)]]$p0),
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
legend_p100 <- c(paste("n =", floor(p * ratios)), expression(paste("true sparsity / default ", p[0])))
x_label <- expression(paste("Prior guess of the number of connections, ", p[0]))

# Create and save plot (Fig. xx in main manuscript).
#plot_sparsity_mcc_f1(p, mean_scores_p100_sf, true_sparsity = true_sparsity_sf_p100, x_label = x_label,
#                     legend = legend_p100, ylim_sparsity = c(min(est_sparsities_sf_p100), max(est_sparsities_sf_p100)),
#                     filename = "Figures/Supplementary/Tau_analysis_p100_sf.pdf")


plot_sparsity_mcc_sl(p, mean_scores_p100_sf, true_sparsity = true_sparsity_sf_p100, x_label = x_label,
                     legend = legend_p100, ylim_sparsity = c(min(est_sparsities_sf_p100), max(est_sparsities_sf_p100)),
                     ylim_steinloss = c(min(mean_scores_p100_sf$n400$sl_omega), max(mean_scores_p100_sf$n60$sl_omega)),
                     filename = "Figures/Main_article/Tau_analysis_p100_sf.eps", save_eps = TRUE, save_pdf = FALSE,
                     image_size = c(9,3))

# Corresponding tau^2 values.
min_p0 <- min(mean_scores_p100_sf$n60$p0)
max_p0 <- max(mean_scores_p100_sf$n60$p0)
round(cbind(calculate_tau(p * ratios, p, min_p0),
            calculate_tau(p * ratios, p, max_p0)), 4)


### Analysis with p = 100 and random structure with network size p/2.
# Define needed parameters.
p <- 100
ratios <- c(0.6, 0.8, 1, 1.25, 1.5, 2, 4)
n_repeats <- 20
n_p0s <- 75
step_p0 <- 5
structure = "random"
seed <- 22102027

# Run analysis if results file does not exist.
if (!file.exists("Results_files/Tau_analyses/scores_p100_random.rds")) {
  print(s2 <- Sys.time())
  scores_p100_random <- list()
  for (ratio in ratios) {
    n <- floor(p * ratio)
    cat("n =", n, "\n")
    scores_p100_random[[paste0("n", n)]] <- calculation_over_p0(p = p, n = n,
                                                                n_repeats = n_repeats,
                                                                n_p0s = n_p0s,
                                                                step_p0 = step_p0,
                                                                structure = structure,
                                                                seed = seed)
  }
  print(Sys.time() - s2)
  saveRDS(scores_p100_random, file = "Results_files/Tau_analyses/scores_p100_random.rds")
}

# Load results if analysis is already done.
if (file.exists("Results_files/Tau_analyses/scores_p100_random.rds")) {
  scores_p100_random <- readRDS("Results_files/Tau_analyses/scores_p100_random.rds")
}

# Calculate mean scores.
mean_scores_p100_random <- list()
for (ratio in ratios) {
  n <- floor(p * ratio)
  mean_scores_p100_random[[paste0("n", n)]] <- aggregate(scores_p100_random[[paste0("n", n)]],
                                                         by = list(scores_p100_random[[paste0("n", n)]]$p0),
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
legend_p100 <- c(paste("n =", floor(p * ratios)), expression(paste("true sparsity / default ", p[0]))) 
x_label <- expression(paste("Prior guess of the number of connections, ", p[0]))

# Create and save plot.
#plot_sparsity_mcc_f1(p, mean_scores_p100_random, true_sparsity = true_sparsity_random_p100, x_label = x_label,
#                     legend = legend_p100, ylim_sparsity = c(min(est_sparsities_random_p100), max(est_sparsities_random_p100)),
#                     filename = "Figures/Supplementary/Tau_analysis_p100_random.pdf")

# Corresponding tau^2 values.
min_p0 <- min(mean_scores_p100_random$n60$p0)
max_p0 <- max(mean_scores_p100_random$n60$p0)
round(cbind(calculate_tau(p * ratios, p, min_p0),
            calculate_tau(p * ratios, p, max_p0)), 4)


### Analysis with p = 100 and random structure with network size p * 1.5.
# Define needed parameters.
p <- 100
ratios <- c(0.6, 0.8, 1, 1.25, 1.5, 2, 4)
net_size <- floor(p * 1.5)
n_repeats <- 20
n_p0s <- 75
step_p0 <- 5
structure = "random"
seed <- 22102028

# Run analysis if results file does not exist.
if (!file.exists("Results_files/Tau_analyses/scores_p100_random_denser.rds")) {
  print(s3 <- Sys.time())
  scores_p100_random_denser <- list()
  for (ratio in ratios) {
    n <- floor(p * ratio)
    cat("n =", n, "\n")
    scores_p100_random_denser[[paste0("n", n)]] <- calculation_over_p0(p = p, n = n,
                                                                       n_repeats = n_repeats,
                                                                       n_p0s = n_p0s,
                                                                       step_p0 = step_p0,
                                                                       structure = structure,
                                                                       net_size = net_size,
                                                                       seed = seed)
  }
  print(Sys.time() - s3)
  saveRDS(scores_p100_random_denser, file = "Results_files/Tau_analyses/scores_p100_random_denser.rds")
}

# Load results if analysis is already done.
if (file.exists("Results_files/Tau_analyses/scores_p100_random_denser.rds")) {
  scores_p100_random_denser <- readRDS("Results_files/Tau_analyses/scores_p100_random_denser.rds")
}

# Calculate mean scores.
mean_scores_p100_random_denser <- list()
for (ratio in ratios) {
  n <- floor(p * ratio)
  mean_scores_p100_random_denser[[paste0("n", n)]] <- aggregate(scores_p100_random_denser[[paste0("n", n)]],
                                                                by = list(scores_p100_random_denser[[paste0("n", n)]]$p0),
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
legend_p100 <- c(paste("n =", floor(p * ratios)), expression(paste("true sparsity / default ", p[0])))
x_label <- expression(paste("Prior guess of the number of connections, ", p[0]))

# Create and save plot.
#plot_sparsity_mcc_f1(p, mean_scores_p100_random_denser, true_sparsity = true_sparsity_rnd_dens_p100,
#                     x_label = x_label, legend = legend_p100,
#                     ylim_sparsity = c(min(est_sparsities_rnd_dens_p100), max(est_sparsities_rnd_dens_p100)),
#                     filename = "Figures/Supplementary/Tau_analysis_p100_random_denser.pdf")

# Corresponding tau^2 values.
min_p0 <- min(mean_scores_p100_random_denser$n60$p0)
max_p0 <- max(mean_scores_p100_random_denser$n60$p0)
round(cbind(calculate_tau(p * ratios, p, min_p0),
            calculate_tau(p * ratios, p, max_p0)), 4)

### Analysis with p = 100 and hub structure.
# Define needed parameters.
p <- 100
ratios <- c(0.6, 0.8, 1, 1.25, 1.5, 2, 4)
n_repeats <- 20
n_p0s <- 75
step_p0 <- 5
structure = "hub"
seed <- 11122025

# Run analysis if results file does not exist.
if (!file.exists("Results_files/Tau_analyses/scores_p100_hub.rds")) {
  print(s6 <- Sys.time())
  scores_p100_hub <- list()
  for (ratio in ratios) {
    n <- floor(p * ratio)
    cat("n =", n, "\n")
    scores_p100_hub[[paste0("n", n)]] <- calculation_over_p0(p = p, n = n,
                                                             n_repeats = n_repeats,
                                                             n_p0s = n_p0s,
                                                             step_p0 = step_p0,
                                                             structure = structure,
                                                             seed = seed)
  }
  print(Sys.time() - s6)
  saveRDS(scores_p100_hub, file = "Results_files/Tau_analyses/scores_p100_hub.rds")
}

# Load results if analysis is already done.
if (file.exists("Results_files/Tau_analyses/scores_p100_hub.rds")) {
  scores_p100_hub <- readRDS("Results_files/Tau_analyses/scores_p100_hub.rds")
}

# Calculate mean scores.
mean_scores_p100_hub <- list()
for (ratio in ratios) {
  n <- floor(p * ratio)
  mean_scores_p100_hub[[paste0("n", n)]] <- aggregate(scores_p100_hub[[paste0("n", n)]],
                                                      by = list(scores_p100_hub[[paste0("n", n)]]$p0),
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
legend_p100 <- c(paste("n =", floor(p * ratios)), expression(paste("true sparsity / default ", p[0])))
x_label <- expression(paste("Prior guess of the number of connections, ", p[0]))

# Create and save plot.
#plot_sparsity_mcc_f1(p, mean_scores_p100_hub, true_sparsity = true_sparsity_rnd_dens_p100,
#                     x_label = x_label, legend = legend_p100,
#                     ylim_sparsity = c(min(est_sparsities_hub_p100), max(est_sparsities_hub_p100)),
#                     filename = "Figures/Supplementary/Tau_analysis_p100_hub.pdf")

# Corresponding tau^2 values.
min_p0 <- min(mean_scores_p100_hub$n60$p0)
max_p0 <- max(mean_scores_p100_hub$n60$p0)
round(cbind(calculate_tau(p * ratios, p, min_p0),
            calculate_tau(p * ratios, p, max_p0)), 4)

### Combine all plots with p = 100 in one figure.
# pdf("Figures/Supplementary/Tau_analysis_p100.pdf", width = 10.5, height = 10.5)
# par(mfrow = c(3,3))
# plot_sparsity_mcc_f1(100, mean_scores_p100_sf, true_sparsity = true_sparsity_sf_p100,
#                      x_label = x_label, legend = legend_p100, plot_rows = 3,
#                      ylim_sparsity = c(min(est_sparsities_sf_p100), max(est_sparsities_sf_p100)))
# plot_sparsity_mcc_f1(100, mean_scores_p100_random, true_sparsity = true_sparsity_random_p100,
#                      x_label = x_label, legend = legend_p100, plot_rows = 3,
#                      ylim_sparsity = c(min(est_sparsities_random_p100), max(est_sparsities_random_p100)))
# plot_sparsity_mcc_f1(100, mean_scores_p100_random_denser, true_sparsity = true_sparsity_rnd_dens_p100,
#                      x_label = x_label, legend = legend_p100, plot_rows = 3,
#                      ylim_sparsity = c(min(est_sparsities_rnd_dens_p100), max(est_sparsities_rnd_dens_p100)))
# 
# dev.off()

### Combine all plots with p = 100 in one figure.
pdf("Figures/Supplementary/Tau_analysis_with_sl_p100.pdf", width = 10.5, height = 10.5)
par(mfrow = c(3,3))
plot_sparsity_mcc_sl(100, mean_scores_p100_random, true_sparsity = true_sparsity_random_p100,
                     x_label = x_label, legend = legend_p100, plot_rows = 3,
                     ylim_sparsity = c(min(est_sparsities_random_p100), max(est_sparsities_random_p100)),
                     ylim_steinloss = c(min(mean_scores_p100_random$n400$sl_omega), max(mean_scores_p100_random$n60$sl_omega)))
plot_sparsity_mcc_sl(100, mean_scores_p100_random_denser, true_sparsity = true_sparsity_rnd_dens_p100,
                     x_label = x_label, legend = legend_p100, plot_rows = 3,
                     ylim_sparsity = c(min(est_sparsities_rnd_dens_p100), max(est_sparsities_rnd_dens_p100)),
                     ylim_steinloss = c(min(mean_scores_p100_random_denser$n400$sl_omega), max(mean_scores_p100_random_denser$n60$sl_omega)))
plot_sparsity_mcc_sl(100, mean_scores_p100_hub, true_sparsity = true_sparsity_hub_p100,
                     x_label = x_label, legend = legend_p100, plot_rows = 3,
                     ylim_sparsity = c(min(est_sparsities_hub_p100), max(est_sparsities_hub_p100)),
                     ylim_steinloss = c(min(mean_scores_p100_hub$n400$sl_omega), max(mean_scores_p100_hub$n60$sl_omega)))

dev.off()

mean_scores_p100_hub

### Analysis with p = 200 and hub structure.
# Define needed parameters.
p <- 200
ratios <- c(0.6, 0.8, 1, 1.25, 1.5, 2, 4)
n_repeats <- 10
p0s <- c(seq(49, 240, 10))
structure = "hub"
seed <- 22102026

# Run analysis if results file does not exist.
if (!file.exists("Results_files/Tau_analyses/scores_p200_hub.rds")) {
  print(s4 <- Sys.time())
  scores_p200_hub <- list()
  for (ratio in ratios) {
    n <- floor(p * ratio)
    cat("n =", n, "\n")
    scores_p200_hub[[paste0("n", n)]] <- calculation_over_p0(p = p, n = n,
                                                             n_repeats = n_repeats,
                                                             p0s = p0s,
                                                             structure = structure,
                                                             seed = seed)
  }
  print(Sys.time() - s4)
  saveRDS(scores_p200_hub, file = "Results_files/Tau_analyses/scores_p200_hub.rds")
}

# Load results if analysis is already done.
if (file.exists("Results_files/Tau_analyses/scores_p200_hub.rds")) {
  scores_p200_hub <- readRDS("Results_files/Tau_analyses/scores_p200_hub.rds")
}

# Calculate mean scores.
mean_scores_p200_hub <- list()
for (ratio in ratios) {
  n <- floor(p * ratio)
  mean_scores_p200_hub[[paste0("n", n)]] <- aggregate(scores_p200_hub[[paste0("n", n)]],
                                                      by = list(scores_p200_hub[[paste0("n", n)]]$p0),
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
legend_p200 <- c(paste("n =", floor(p * ratios)), expression(paste("true sparsity / default ", p[0]))) 
x_label <- expression(paste("Prior guess of the number of connections, ", p[0]))

# # Create and save plot.
# plot_sparsity_mcc_f1(p, mean_scores_p200_hub, true_sparsity = true_sparsity_hub_p200, x_label = x_label,
#                      legend = legend_p200, legend_position = "topright",
#                      ylim_sparsity = c(min(est_sparsities_hub_p200), max(est_sparsities_hub_p200)),
#                      filename = "Figures/Supplementary/Tau_analysis_p200_hub.pdf")

# Corresponding tau^2 values.
min_p0 <- min(mean_scores_p200_hub$n120$p0)
max_p0 <- max(mean_scores_p200_hub$n120$p0)
round(cbind(calculate_tau(p * ratios, p, min_p0),
            calculate_tau(p * ratios, p, max_p0)), 5)


### Analysis with p = 250 and cluster structure.
# Define needed parameters.
p <- 250
ratios <- c(0.6, 0.8, 1, 1.25, 1.5, 2, 4)
n_repeats <- 10
p0s <- seq(p - 1 - 25*4, p - 1 + 25 * 12, 25)
structure = "cluster"
seed <- 22102029

# Run analysis if results file does not exist.
if (!file.exists("Results_files/Tau_analyses/scores_p250_cluster.rds")) {
  print(s5 <- Sys.time())
  scores_p250_cluster <- list()
  for (ratio in ratios) {
    n <- floor(p * ratio)
    cat("n =", n, "\n")
    scores_p250_cluster[[paste0("n", n)]] <- calculation_over_p0(p = p, n = n,
                                                                 n_repeats = n_repeats,
                                                                 p0s = p0s,
                                                                 structure = structure,
                                                                 net_size = net_size,
                                                                 seed = seed)
  }
  print(Sys.time() - s5)
  saveRDS(scores_p250_cluster, file = "Results_files/Tau_analyses/scores_p250_cluster.rds")
}

# Load results if analysis is already done.
if (file.exists("Results_files/Tau_analyses/scores_p250_cluster.rds")) {
  scores_p250_cluster <- readRDS("Results_files/Tau_analyses/scores_p250_cluster.rds")
}

# Calculate mean scores.
mean_scores_p250_cluster <- list()
for (ratio in ratios) {
  n <- floor(p * ratio)
  mean_scores_p250_cluster[[paste0("n", n)]] <- aggregate(scores_p250_cluster[[paste0("n", n)]],
                                                          by = list(scores_p250_cluster[[paste0("n", n)]]$p0),
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
legend_p250 <- c(paste("n =", floor(p * ratios)), expression(paste("true sparsity / default ", p[0]))) 
x_label <- expression(paste("Prior guess of the number of connections, ", p[0]))

# Create and save plot.
plot_sparsity_mcc_f1(p, mean_scores_p250_cluster, true_sparsity = true_sparsity_cluster_p250, x_label = x_label,
                     legend = legend, legend_position = "bottomright", legend_in_plot = 2,
                     ylim_sparsity = c(min(est_sparsities_cluster_p250), max(est_sparsities_cluster_p250)),
                     filename = "Figures/Supplementary/Tau_analysis_p250_cluster.pdf")

# Corresponding tau^2 values.
min_p0 <- min(mean_scores_p250_cluster$n150$p0)
max_p0 <- max(mean_scores_p250_cluster$n150$p0)
round(cbind(calculate_tau(p * ratios, p, min_p0),
            calculate_tau(p * ratios, p, max_p0)), 5)


### Combine plots with p = 200 and 250 in one figure.
# pdf("Figures/Supplementary/Tau_analysis_p200-250.pdf", width = 10.5, height = 7)
# par(mfrow = c(2,3))
# plot_sparsity_mcc_f1(200, mean_scores_p200_hub, true_sparsity = true_sparsity_hub_p200,
#                      x_label = x_label, legend = legend_p200, plot_rows = 2, legend_position = "topright",
#                      ylim_sparsity = c(min(est_sparsities_hub_p200), max(est_sparsities_hub_p200)))
# plot_sparsity_mcc_f1(250, mean_scores_p250_cluster, true_sparsity = true_sparsity_cluster_p250,
#                      x_label = x_label, legend = legend_p250, plot_rows = 2,
#                      legend_in_plot = 2, legend_position = "bottomright",
#                      ylim_sparsity = c(min(est_sparsities_cluster_p250), max(est_sparsities_cluster_p250)))
# dev.off()

### Combine plots with p = 200 and 250 in one figure.
pdf("Figures/Supplementary/Tau_analysis_with_sl_p200-250.pdf", width = 10.5, height = 7)
par(mfrow = c(2,3))
plot_sparsity_mcc_sl(200, mean_scores_p200_hub, true_sparsity = true_sparsity_hub_p200,
                     x_label = x_label, legend = legend_p200, plot_rows = 2, legend_position = "topright",
                     ylim_sparsity = c(min(est_sparsities_hub_p200), max(est_sparsities_hub_p200)),
                     ylim_steinloss = c(min(mean_scores_p200_hub$n800$sl_omega), max(mean_scores_p200_hub$n120$sl_omega)))
plot_sparsity_mcc_sl(250, mean_scores_p250_cluster, true_sparsity = true_sparsity_cluster_p250,
                     x_label = x_label, legend = legend_p250, plot_rows = 2,
                     legend_in_plot = 2, legend_position = "bottomright",
                     ylim_sparsity = c(min(est_sparsities_cluster_p250), max(est_sparsities_cluster_p250)),
                     ylim_steinloss = c(min(mean_scores_p250_cluster$n1000$sl_omega), max(mean_scores_p250_cluster$n150$sl_omega)))
dev.off()
