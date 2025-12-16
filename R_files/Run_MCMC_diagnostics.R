library(posterior)
library(R.matlab)
library(doParallel)

create_boxplots <- function(post_summaries, p, score, warm_ups, mcmc_lengths, ylim = c(0,1),
                            label_position = -0.03, plot_one_row = TRUE, save_pdf = FALSE,
                            save_png = FALSE, filename = NULL, image_size = c(10, 7), res = 600) {
  if (!is.null(filename) & save_pdf) pdf(filename, width = image_size[1], height = image_size[2])
  else if (!is.null(filename) & save_png) png(filename, width = image_size[1]*res, height = image_size[2]*res,
                                              res = res)
  scores_per_warmup <- list()
  for (wup in warm_ups) {
    scores_per_warmup[[paste0(wup)]] <- matrix(nrow = nrow(post_summaries[[paste0(wup)]][[1]]),
                                               ncol = length(post_summaries[[paste0(wup)]]))
    for (i in 1:length(post_summaries[[paste0(wup)]])) {
      scores_per_warmup[[paste0(wup)]][,i] <- as.matrix(post_summaries[[paste0(wup)]][[i]][, score])
    }
    colnames(scores_per_warmup[[paste0(wup)]]) <- mcmc_lengths
  }
  if(plot_one_row) par(mfrow = c(1,3))
  par(mar = c(3.1, 3, 1.2, 0.2))
  par(mgp = c(2,1,0))
  for (wup in warm_ups) {
    if (wup == 0) title <- "No warmup"
    else if (wup == 500) title <- "500 MCMC samples discarded at the start"
    else if (wup == 1000) title <- "1000 MCMC samples discarded at the start"
    if (score == "rhat") xlabel <- "Rhat"
    else if (score == "ess_bulk") xlabel <- "ESS"
    boxplot(scores_per_warmup[[paste0(wup)]], xaxt = "n", yaxt = "n", xlab = xlabel,
            ylim = ylim, main = title, horizontal = TRUE, at = rev(1:length(mcmc_lengths)))
    
    axis(side = 2, las = 2, labels = FALSE, at = 1:length(mcmc_lengths))
    axis(side = 1)
    
    text(y = 1:length(mcmc_lengths),
         x = par("usr")[1] + label_position * (par("usr")[2] - par("usr")[1]),
         labels = rev(mcmc_lengths),
         xpd = NA,
         cex = 1, adj = 1)
    grid(ny = NA)
  }
  if (!is.null(filename) & (save_pdf | save_png)) dev.off()
}

read_matlab_matrices <- function(path, datanros, n_mcmc) {
  all_matrices <- list()
  for (datanro in datanros) {
    cat("Reading dataset ", datanro, "...", sep = "")
    matrices <- readMat(paste0(path_matrices, structure, "_p", p, "_Omega_samples_datanro_", datanro, ".mat"))
    if (p == 100) {
      ltri_matrices <- matrix(nrow = n_mcmc, ncol = p*(p-1)/2+p)
      for (i in 1:n_mcmc) {
        ltri_matrices[i,] <- matrices$Omega[,,i][lower.tri(diag(p), diag = TRUE)]
      }
      all_matrices[[paste0(datanro)]] <- ltri_matrices
    }
    else all_matrices[[paste0(datanro)]] <- matrices$Omega
    cat("done.\n")
  }
  gc()
  return(all_matrices)
}

calculate_posterior_summaries <- function(all_matrices, datanros, n_warmups = c(0, 500, 1000),
                                          n_mcmc_samples = c(1000, 2000, 3000, 4000, 5000, 6000),
                                          n_cores = 0) {
  if (n_cores == 0) n_cores <- min(c(length(n_mcmc_samples)*length(n_warmups), detectCores() - 1))
  posterior_summaries <- list()
  for (datanro in datanros) {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    summaries_per_data <- foreach (n_warmup = n_warmups) %:%
      foreach (i = 1:length(n_mcmc_samples), .packages = c("posterior")) %dopar% {
        n_mcmc_sample <- n_mcmc_samples[i]
        summarise_draws(all_matrices[[paste0(datanro)]][(1+n_warmup):(n_warmup + n_mcmc_sample),])
      }
    stopCluster(cl)
    stopImplicitCluster()
    posterior_summaries[[paste0(datanro)]] <- summaries_per_data
    names(posterior_summaries[[paste0(datanro)]]) <- as.character(n_warmups)
    for (n_warmup in n_warmups) names(posterior_summaries[[paste0(datanro)]][[paste0(n_warmup)]]) <- as.character(n_mcmc_samples)
  }
  return(posterior_summaries)
}

combine_over_datasets <- function(posterior_summaries, datanros, n_warmups, n_mcmc_samples) {
  combined_scores <- list()
  for (n_warmup in n_warmups) {
    for (n_mcmc in n_mcmc_samples) {
      for (datanro in datanros) {
        combined_scores[[paste0(n_warmup)]][[paste0(n_mcmc)]] <- rbind(combined_scores[[paste0(n_warmup)]][[paste0(n_mcmc)]],
                                                                       posterior_summaries[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(n_mcmc)]])
      }
    }
  }
  return(combined_scores)
}

convert_array_inds_to_ltri_ind <- function(indices) {
  indices <- sort(indices, decreasing = TRUE)
  A <- matrix(1:p^2, ncol = p, nrow = p, byrow = TRUE)
  ind <- which(A[lower.tri(A, diag = TRUE)] == A[indices[1], indices[2]])
  return(ind)
}

create_adj_mat_CI <- function(samples, p, quant = 0.25) {
  low_q <- up_q <- matrix(0, p , p)
  low_q[lower.tri(diag(p), diag = TRUE)] <- apply(samples, 2, quantile, quant)
  up_q[lower.tri(diag(p), diag = TRUE)] <- apply(samples, 2, quantile, 1 - quant)
  Theta_est <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:i) {
      if (i == j) {
        next
      }
      else if (low_q[i,j] <= 0 & up_q[i,j] > 0) {
        next
      }
      else {
        Theta_est[i,j] <- 1
        Theta_est[j,i] <- 1
      }
    }
  }
  return(Theta_est)
}

#####
### Datasets with p = 100.

### MCMC diagnostics for random network datasets generated using bdgraph with p = 100.
# Initialize parameters.
p <- 100
structure <- "random"
path_matrices <- paste0("Results_files/p", p, "/GHS_MCMC/",  structure,  "/")
datanros <- c(9, 25, 33, 42, 47)
n_warmups <- c(0, 500, 1000)
mcmc_lengths <- seq(1000, 6000, 1000)

# Read MCMC samples from the MATLAB files.
random_5_datasets <- read_matlab_matrices(path, datanros, n_mcmc = 10000)

# Calculate posterior summaries.
random_post_summaries <- calculate_posterior_summaries(random_5_datasets, datanros = datanros,
                                                       n_cores = 16)
# Combine scores over datasets.
random_combined_scores <- combine_over_datasets(random_post_summaries, datanros, n_warmups, mcmc_lengths)

# Save only rhat, ess_bulk and ess_tail.
for (datanro in datanros) {
  for (n_warmup in n_warmups) {
    for (mcmc_length in mcmc_lengths) {
      random_post_summaries[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(mcmc_length)]] <- random_post_summaries[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(mcmc_length)]][,c("rhat", "ess_bulk", "ess_tail")]
    }
  }
}
for (n_warmup in n_warmups) {
  for (mcmc_length in mcmc_lengths) {
    random_combined_scores[[paste0(n_warmup)]][[paste0(mcmc_length)]] <- random_combined_scores[[paste0(n_warmup)]][[paste0(mcmc_length)]][,c("rhat", "ess_bulk", "ess_tail")]
  }
}

# Save summaries.
saveRDS(random_post_summaries, file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/random_per_dataset.Rds"))
saveRDS(random_combined_scores, file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/random_combined.Rds"))

# Read summaries if already done.
random_post_summaries <- readRDS(file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/random_per_dataset.Rds"))
random_combined_scores <- readRDS(file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/random_combined.Rds"))

# Plot per dataset boxplots for Rhats and ESSs.
par(mfrow = c(5,3))
for (datanro in datanros) {
  create_boxplots(random_post_summaries[[paste0(datanro)]], p, "rhat", warm_ups = n_warmups,
                  mcmc_lengths = mcmc_lengths, ylim = c(1, 1.2), plot_one_row = FALSE)
}
par(mfrow = c(5,3))
for (datanro in datanros) {
  create_boxplots(random_post_summaries[[paste0(datanro)]], p, "ess_bulk", warm_ups = n_warmups,
                  mcmc_lengths = mcmc_lengths, ylim = c(0, 7000), plot_one_row = FALSE)
}

# Create boxplots for combined results.
create_boxplots(random_combined_scores, p, "rhat", n_warmups, mcmc_lengths, ylim = c(1, 1.2),
                filename = "Figures/Supplementary/random_p100_Rhats.png", image_size = c(11, 5),
                save_png = TRUE, res = 800)

create_boxplots(random_combined_scores, p, "ess_bulk", n_warmups, mcmc_lengths, ylim = c(0, 7000),
                filename = "Figures/Supplementary/random_p100_ESSs.png", image_size = c(11, 5),
                save_png = TRUE, res = 800)


### MCMC diagnostics for scale-free network datasets generated using bdgraph with p = 100.
# Initialize parameters.
p <- 100
structure <- "bdgraph_sf"
path_matrices <- paste0("Results_files/p", p, "/GHS_MCMC/",  structure,  "/")
datanros <- c(9, 25, 33, 42, 47)
n_warmups <- c(0, 500, 1000)
mcmc_lengths <- seq(1000, 6000, 1000)

# For trace plots, read only dataset 25.
bdgraph_sf_dataset_25_p100 <- read_matlab_matrices(path, 25, n_mcmc = 10000)

# Load simulations to check zeros and non-zeros.
path_data <- paste0("Data/n", 120, "_p", p, "/")
sim_bdgraph_sf_p100 <- readRDS(file = paste0(path_data, "bdgraph_scale-free_n", 120, "_p", p, ".Rds"))

sim <- sim_bdgraph_sf_p100[[25]]
which(sim$theta != 0, arr.ind = TRUE)

diagonals <- list(c(1,1), c(2,2), c(5,5), c(31,31), c(43,43), c(76,76))
non_zeros <- list(c(1,2), c(2,24), c(5,95), c(35,30), c(53,10), c(73,65))
zeros <- list(c(1,10), c(2,9), c(5,23), c(35,65), c(53,55), c(72,98))

# Start of the trace plots
pdf("Figures/Supplementary/Bdgraph_sf_data_25_traceplots.pdf", width = 7, 8.5)
cols <- palette.colors(6)
ltys <- rep(c(2,4,5,6,1), length.out = 6)
par(mfrow = c(3,1), mar = c(3.1, 3.1, 0.2, 0.2), mgp = c(2,1,0))
plot(bdgraph_sf_dataset_25_p100[["25"]][1:1500, convert_array_inds_to_ltri_ind(diagonals[[1]])], type = "l",
     ylim = c(0, 27), lty = ltys[1], xlab = "Iteration", ylab = "Diagonal of the precision matrix")
for (i in 2:length(diagonals)) {
  lines(bdgraph_sf_dataset_25_p100[["25"]][1:1500, convert_array_inds_to_ltri_ind(diagonals[[i]])],
        col = cols[i], lty = ltys[i])
}

plot(bdgraph_sf_dataset_25_p100[["25"]][1:1500, convert_array_inds_to_ltri_ind(non_zeros[[1]])], type = "l",
     ylim = c(-9, 0.5), lty = ltys[1], xlab = "Iteration",
     ylab = "Non-zero off-diagonal of the precision matrix")
for (i in 2:length(non_zeros)) {
  lines(bdgraph_sf_dataset_25_p100[["25"]][1:1500, convert_array_inds_to_ltri_ind(non_zeros[[i]])],
        col = cols[i], lty = ltys[i])
}

plot(bdgraph_sf_dataset_25_p100[["25"]][1:1500, convert_array_inds_to_ltri_ind(zeros[[1]])], type = "l",
     ylim = c(-1.5, 1), lty = ltys[1], xlab = "Iteration",
     ylab = "Zero off-diagonal of the precision matrix")
for (i in 2:length(non_zeros)) {
  lines(bdgraph_sf_dataset_25_p100[["25"]][1:1500,convert_array_inds_to_ltri_ind(zeros[[i]])],
        col = cols[i], lty = ltys[i])
}
dev.off()

# End of the trace plots
pdf("Figures/Supplementary/Bdgraph_sf_data_42_traceplots_end.pdf", width = 7, 8.5)
cols <- palette.colors(6)
ltys <- rep(c(2,4,5,6,1), length.out = 6)
par(mfrow = c(3,1), mar = c(3.1, 3.1, 0.2, 0.2), mgp = c(2,1,0))
plot(1501:10000, bdgraph_sf_dataset_25_p100[["25"]][1501:10000, convert_array_inds_to_ltri_ind(diagonals[[1]])],
     type = "l", ylim = c(0, 27), lty = ltys[1], xlab = "Iteration", ylab = "Diagonal of the precision matrix")
for (i in 2:length(diagonals)) {
  lines(1501:10000, bdgraph_sf_dataset_25_p100[["25"]][1501:10000, convert_array_inds_to_ltri_ind(diagonals[[i]])],
        col = cols[i], lty = ltys[i])
}
plot(1501:10000, bdgraph_sf_dataset_25_p100[["25"]][1501:10000, convert_array_inds_to_ltri_ind(non_zeros[[1]])],
     type = "l", ylim = c(-9, 0.5), lty = ltys[1], xlab = "Iteration",
     ylab = "Non-zero off-diagonal of the precision matrix")
for (i in 2:length(non_zeros)) {
  lines(1501:10000, bdgraph_sf_dataset_25_p100[["25"]][1501:10000, convert_array_inds_to_ltri_ind(non_zeros[[i]])],
        col = cols[i], lty = ltys[i])
}
plot(1501:10000, bdgraph_sf_dataset_25_p100[["25"]][1501:10000, convert_array_inds_to_ltri_ind(zeros[[1]])],
     type = "l", ylim = c(-1.5, 1), lty = ltys[1], xlab = "Iteration",
     ylab = "Zero off-diagonal of the precision matrix")
for (i in 2:length(non_zeros)) {
  lines(1501:10000, bdgraph_sf_dataset_25_p100[["25"]][1501:10000, convert_array_inds_to_ltri_ind(zeros[[i]])],
        col = cols[i], lty = ltys[i])
}
dev.off()

# Remove datasets to save memory for further analysis.
rm(bdgraph_sf_dataset_25_p100)
rm(sim_bdgraph_sf_p100)
gc()

# Read MCMC samples from the MATLAB files.
bdgraph_sf_5_datasets <- read_matlab_matrices(path, datanros, n_mcmc = 10000)

# Calculate posterior summaries.
bdgraph_sf_post_summaries <- calculate_posterior_summaries(bdgraph_sf_5_datasets, datanros = datanros,
                                                           n_cores = 16)
# Combine scores over datasets.
bdgraph_sf_combined_scores <- combine_over_datasets(bdgraph_sf_post_summaries, "rhat", datanros, n_warmups, mcmc_lengths)

# Save only rhat, ess_bulk and ess_tail.
for (datanro in datanros) {
  for (n_warmup in n_warmups) {
    for (mcmc_length in mcmc_lengths) {
      bdgraph_sf_post_summaries[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(mcmc_length)]] <- bdgraph_sf_post_summaries[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(mcmc_length)]][,c("rhat", "ess_bulk", "ess_tail")]
    }
  }
}
for (n_warmup in n_warmups) {
  for (mcmc_length in mcmc_lengths) {
    bdgraph_sf_combined_scores[[paste0(n_warmup)]][[paste0(mcmc_length)]] <- bdgraph_sf_combined_scores[[paste0(n_warmup)]][[paste0(mcmc_length)]][,c("rhat", "ess_bulk", "ess_tail")]
  }
}

# Save summaries.
saveRDS(bdgraph_sf_post_summaries, file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/bdgraph_sf_per_dataset.Rds"))
saveRDS(bdgraph_sf_combined_scores, file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/bdgraph_sf_combined.Rds"))

# Read summaries if already done.
bdgraph_sf_post_summaries <- readRDS(file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/bdgraph_sf_per_dataset.Rds"))
bdgraph_sf_combined_scores <- readRDS(file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/bdgraph_sf_combined.Rds"))

# Plot per dataset boxplots for Rhats and ESSs.
par(mfrow = c(5,3))
for (datanro in datanros) {
  create_boxplots(bdgraph_sf_post_summaries[[paste0(datanro)]], p, "rhat", warm_ups = n_warmups,
                  mcmc_lengths = mcmc_lengths, ylim = c(1, 1.2), plot_one_row = FALSE)
}
par(mfrow = c(5,3))
for (datanro in datanros) {
  create_boxplots(bdgraph_sf_post_summaries[[paste0(datanro)]], p, "ess_bulk", warm_ups = n_warmups,
                  mcmc_lengths = mcmc_lengths, ylim = c(0, 7000), plot_one_row = FALSE)
}

# Create boxplots for combined results.
create_boxplots(bdgraph_sf_combined_scores, p, "rhat", n_warmups, mcmc_lengths, ylim = c(1, 1.2),
                filename = "Figures/Supplementary/bdgraph_sf_p100_Rhats.png", image_size = c(11, 5),
                save_png = TRUE, res = 800)

create_boxplots(bdgraph_sf_combined_scores, p, "ess_bulk", n_warmups, mcmc_lengths, ylim = c(0, 7000),
                filename = "Figures/Supplementary/bdgraph_sf_p100_ESSs.png", image_size = c(11, 5),
                save_png = TRUE, res = 800)


### MCMC diagnostics for scale-free datasets generated using huge with p = 100.
# Initialize parameters.
p <- 100
structure <- "huge_sf"
path_matrices <- paste0("Results_files/p", p, "/GHS_MCMC/",  structure,  "/")
datanros <- c(9, 25, 33, 42, 47)
n_warmups <- c(0, 500, 1000)
mcmc_lengths <- seq(1000, 6000, 1000)

# Read MCMC samples from the MATLAB files.
huge_sf_5_datasets <- read_matlab_matrices(path, datanros, n_mcmc = 10000)

# Calculate posterior summaries.
huge_sf_post_summaries <- calculate_posterior_summaries(huge_sf_5_datasets, datanros = datanros,
                                                        n_cores = 16)
# Combine scores over datasets.
huge_sf_combined_scores <- combine_over_datasets(huge_sf_post_summaries, "rhat", datanros, n_warmups, mcmc_lengths)

# Save only rhat, ess_bulk and ess_tail.
for (datanro in datanros) {
  for (n_warmup in n_warmups) {
    for (mcmc_length in mcmc_lengths) {
      huge_sf_post_summaries[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(mcmc_length)]] <- huge_sf_post_summaries[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(mcmc_length)]][,c("rhat", "ess_bulk", "ess_tail")]
    }
  }
}
for (n_warmup in n_warmups) {
  for (mcmc_length in mcmc_lengths) {
    huge_sf_combined_scores[[paste0(n_warmup)]][[paste0(mcmc_length)]] <- huge_sf_combined_scores[[paste0(n_warmup)]][[paste0(mcmc_length)]][,c("rhat", "ess_bulk", "ess_tail")]
  }
}

# Save summaries.
saveRDS(huge_sf_post_summaries, file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/huge_sf_per_dataset.Rds"))
saveRDS(huge_sf_combined_scores, file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/huge_sf_combined.Rds"))

# Read summaries if already done.
huge_sf_post_summaries <- readRDS(file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/huge_sf_per_dataset.Rds"))
huge_sf_combined_scores <- readRDS(file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/huge_sf_combined.Rds"))

# Plot per dataset boxplots for Rhats and ESSs.
par(mfrow = c(5,3))
for (datanro in datanros) {
  create_boxplots(huge_sf_post_summaries[[paste0(datanro)]], p, "rhat", warm_ups = n_warmups,
                  mcmc_lengths = mcmc_lengths, ylim = c(1, 1.2), plot_one_row = FALSE)
}
par(mfrow = c(5,3))
for (datanro in datanros) {
  create_boxplots(huge_sf_post_summaries[[paste0(datanro)]], p, "ess_bulk", warm_ups = n_warmups,
                  mcmc_lengths = mcmc_lengths, ylim = c(0, 7000), plot_one_row = FALSE)
}

# Create boxplots for combined results.
create_boxplots(combined_scores, p, "rhat", n_warmups, mcmc_lengths, ylim = c(1, 1.2),
                filename = "Figures/Supplementary/huge_sf_p100_Rhats.png", image_size = c(11, 5),
                save_png = TRUE, res = 800)

create_boxplots(combined_scores, p, "ess_bulk", n_warmups, mcmc_lengths, ylim = c(0, 7000),
                filename = "Figures/Supplementary/huge_sf_p100_ESSs.png", image_size = c(11, 5),
                save_png = TRUE, res = 800)


### MCMC diagnostics for hub network datasets generated using huge with p = 100.
# Initialize parameters.
p <- 100
structure <- "hubs"
path_matrices <- paste0("Results_files/p", p, "/GHS_MCMC/",  structure,  "/")
datanros <- c(9, 25, 33, 42, 47)
n_warmups <- c(0, 500, 1000)
mcmc_lengths <- seq(1000, 6000, 1000)

# Read MCMC samples from the MATLAB files.
hub_5_datasets <- read_matlab_matrices(path, datanros, n_mcmc = 10000)

# Calculate posterior summaries.
hub_post_summaries <- calculate_posterior_summaries(hub_5_datasets, datanros = datanros,
                                                    n_cores = 16)
# Combine scores over datasets.
hub_combined_scores <- combine_over_datasets(hub_post_summaries, "rhat", datanros, n_warmups, mcmc_lengths)

# Save only rhat, ess_bulk and ess_tail.
for (datanro in datanros) {
  for (n_warmup in n_warmups) {
    for (mcmc_length in mcmc_lengths) {
      hub_post_summaries[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(mcmc_length)]] <- hub_post_summaries[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(mcmc_length)]][,c("rhat", "ess_bulk", "ess_tail")]
    }
  }
}
for (n_warmup in n_warmups) {
  for (mcmc_length in mcmc_lengths) {
    hub_combined_scores[[paste0(n_warmup)]][[paste0(mcmc_length)]] <- hub_combined_scores[[paste0(n_warmup)]][[paste0(mcmc_length)]][,c("rhat", "ess_bulk", "ess_tail")]
  }
}

# Save summaries.
saveRDS(hub_post_summaries, file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/hub_per_dataset.Rds"))
saveRDS(hub_combined_scores, file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/hub_combined.Rds"))

# Read summaries if already done.
hub_post_summaries <- readRDS(file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/hub_per_dataset.Rds"))
hub_combined_scores <- readRDS(file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/hub_combined.Rds"))

# Plot per dataset boxplots for Rhats and ESSs.
par(mfrow = c(5,3))
for (datanro in datanros) {
  create_boxplots(hub_post_summaries[[paste0(datanro)]], p, "rhat", warm_ups = n_warmups,
                  mcmc_lengths = mcmc_lengths, ylim = c(1, 1.2), plot_one_row = FALSE)
}
par(mfrow = c(5,3))
for (datanro in datanros) {
  create_boxplots(hub_post_summaries[[paste0(datanro)]], p, "ess_bulk", warm_ups = n_warmups,
                  mcmc_lengths = mcmc_lengths, ylim = c(0, 7000), plot_one_row = FALSE)
}

# Create boxplots for combined results.
create_boxplots(hub_combined_scores, p, "rhat", n_warmups, mcmc_lengths, ylim = c(1, 1.2),
                filename = "Figures/Supplementary/hubs_p100_Rhats.png", image_size = c(11, 5),
                save_png = TRUE, res = 800)

create_boxplots(hub_combined_scores, p, "ess_bulk", n_warmups, mcmc_lengths, ylim = c(0, 7000),
                filename = "Figures/Supplementary/hubs_p100_ESSs.png", image_size = c(11, 5),
                save_png = TRUE, res = 800)

# Remove datasets with p = 100.
rm(random_5_datasets)
rm(bdgraph_sf_5_datasets)
rm(huge_sf_5_datasets)
rm(hub_5_datasets)
gc()

#####
### Datasets with p = 200

### MCMC diagnostics for scale-free network datasets generated using bdgraph with p = 100.
# Initialize parameters.
p <- 200
structure <- "random"
path_matrices <- paste0("Results_files/p", p, "/GHS_MCMC/",  structure,  "/")
datanros <- c(9, 25, 33, 42, 47)
n_warmups <- c(0, 500, 1000)
mcmc_lengths <- seq(1000, 6000, 1000)

# Read MCMC samples from the MATLAB files.
random_5_datasets_p200 <- read_matlab_matrices(path, datanros, n_mcmc = 10000)

# Calculate posterior summaries.
random_post_summaries_p200 <- calculate_posterior_summaries(random_5_datasets_p200, datanros = datanros,
                                                            n_cores = 6)
# Combine scores over datasets.
random_combined_scores_p200 <- combine_over_datasets(random_post_summaries_p200, "rhat",
                                                     datanros, n_warmups, mcmc_lengths)

# Save only rhat, ess_bulk and ess_tail.
for (datanro in datanros) {
  for (n_warmup in n_warmups) {
    for (mcmc_length in mcmc_lengths) {
      random_post_summaries_p200[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(mcmc_length)]] <- random_post_summaries_p200[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(mcmc_length)]][,c("rhat", "ess_bulk", "ess_tail")]
    }
  }
}
for (n_warmup in n_warmups) {
  for (mcmc_length in mcmc_lengths) {
    random_combined_scores_p200[[paste0(n_warmup)]][[paste0(mcmc_length)]] <- random_combined_scores_p200[[paste0(n_warmup)]][[paste0(mcmc_length)]][,c("rhat", "ess_bulk", "ess_tail")]
  }
}

# Save summaries.
saveRDS(random_post_summaries_p200, file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/random_per_dataset_p200.Rds"))
saveRDS(random_combined_scores_p200, file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/random_combined_p200.Rds"))

# Read summaries if already done.
random_post_summaries_p200 <- readRDS(file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/random_per_dataset_p200.Rds"))
random_combined_scores_p200 <- readRDS(file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/random_combined_p200.Rds"))


# Plot per dataset boxplots for Rhats and ESSs.
par(mfrow = c(5,3))
for (datanro in datanros) {
  create_boxplots(random_post_summaries_p200[[paste0(datanro)]], p, "rhat", warm_ups = n_warmups,
                  mcmc_lengths = mcmc_lengths, ylim = c(1, 1.2), plot_one_row = FALSE)
}
par(mfrow = c(5,3))
for (datanro in datanros) {
  create_boxplots(random_post_summaries_p200[[paste0(datanro)]], p, "ess_bulk", warm_ups = n_warmups,
                  mcmc_lengths = mcmc_lengths, ylim = c(0, 7000), plot_one_row = FALSE)
}

# Create boxplots for combined results.
create_boxplots(random_combined_scores_p200, p, "rhat", n_warmups, mcmc_lengths, ylim = c(1, 1.2),
                filename = "Figures/Supplementary/random_p200_Rhats.png", image_size = c(11, 5),
                save_png = TRUE, res = 800)

create_boxplots(random_combined_scores_p200, p, "ess_bulk", n_warmups, mcmc_lengths, ylim = c(0, 7000),
                filename = "Figures/Supplementary/random_p200_ESSs.png", image_size = c(11, 5),
                save_png = TRUE, res = 800)

# Remove datasets to save memory for further analysis.
rm(random_5_datasets_p200)
gc()

### Not performed!!! MCMC diagnostics for scale-free network datasets generated using bdgraph with p = 100.
# Initialize parameters.
# p <- 200
# structure <- "bdgraph_sf"
# path_matrices <- paste0("Results_files/p", p, "/GHS_MCMC/",  structure,  "/")
# datanros <- c(9, 25, 33, 42, 47)
# n_warmups <- c(0, 500, 1000)
# mcmc_lengths <- seq(1000, 6000, 1000)
# 
# # Read MCMC samples from the MATLAB files.
# bdgraph_sf_5_datasets_p200 <- read_matlab_matrices(path, datanros, n_mcmc = 10000)
# 
# # Calculate posterior summaries.
# bdgraph_sf_post_summaries_p200 <- calculate_posterior_summaries(bdgraph_sf_5_datasets_p200, datanros = datanros,
#                                                                 n_cores = 16)
# # Combine scores over datasets.
# bdgraph_sf_combined_scores_p200 <- combine_over_datasets(bdgraph_sf_post_summaries_p200, "rhat",
#                                                          datanros, n_warmups, mcmc_lengths)
# 
# Save only rhat, ess_bulk and ess_tail.
# for (datanro in datanros) {
#   for (n_warmup in n_warmups) {
#     for (mcmc_length in mcmc_lengths) {
#       bdgraph_sf_post_summaries_p200[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(mcmc_length)]] <- bdgraph_sf_post_summaries_p200[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(mcmc_length)]][,c("rhat", "ess_bulk", "ess_tail")]
#     }
#   }
# }
# for (n_warmup in n_warmups) {
#   for (mcmc_length in mcmc_lengths) {
#     bdgraph_sf_combined_scores_p200[[paste0(n_warmup)]][[paste0(mcmc_length)]] <- bdgraph_sf_combined_scores_p200[[paste0(n_warmup)]][[paste0(mcmc_length)]][,c("rhat", "ess_bulk", "ess_tail")]
#   }
# }
#
# # Save summaries.
# saveRDS(bdgraph_sf_post_summaries_p200, file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/bdgraph_sf_per_dataset.Rds"))
# saveRDS(bdgraph_sf_combined_scores_p200, file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/bdgraph_sf_combined.Rds"))
# 
# # Read summaries if already done.
# bdgraph_sf_post_summaries_p200 <- readRDS(file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/bdgraph_sf_per_dataset.Rds"))
# bdgraph_sf_combined_scores_p200 <- readRDS(file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/bdgraph_sf_combined.Rds"))
# 
# # Plot per dataset boxplots for Rhats and ESSs.
# par(mfrow = c(5,3))
# for (datanro in datanros) {
#   create_boxplots(bdgraph_sf_post_summaries_p200[[paste0(datanro)]], p, "rhat", warm_ups = n_warmups,
#                   mcmc_lengths = mcmc_lengths, ylim = c(1, 1.2), plot_one_row = FALSE)
# }
# par(mfrow = c(5,3))
# for (datanro in datanros) {
#   create_boxplots(bdgraph_sf_post_summaries_p200[[paste0(datanro)]], p, "ess_bulk", warm_ups = n_warmups,
#                   mcmc_lengths = mcmc_lengths, ylim = c(0, 7000), plot_one_row = FALSE)
# }
# 
# # Create boxplots for combined results.
# create_boxplots(bdgraph_sf_combined_scores_p200, p, "rhat", n_warmups, mcmc_lengths, ylim = c(1, 1.2))
# 
# create_boxplots(bdgraph_sf_combined_scores_p200, p, "ess_bulk", n_warmups, mcmc_lengths, ylim = c(0, 7000))
#
# # Remove datasets to save memory for further analysis.
# rm(bdgraph_sf_5_datasets_p200)
# gc()


### MCMC diagnostics for scale-free network datasets generated using huge with p = 200.
# Initialize parameters.
p <- 200
structure <- "huge_sf"
path_matrices <- paste0("Results_files/p", p, "/GHS_MCMC/",  structure,  "/")
datanros <- c(9, 25, 33, 42, 47)
n_warmups <- c(0, 500, 1000)
mcmc_lengths <- seq(1000, 6000, 1000)

# For trace plots, read only dataset 42.
huge_sf_dataset_42_p200 <- read_matlab_matrices(path, 42, n_mcmc = 10000)

# Load simulations to check zeros and non-zeros.
path_data <- paste0("Data/n", 120, "_p", p, "/")
sim_huge_sf_p200 <- readRDS(file = paste0(path_data, "huge_scale-free_n", 120, "_p", p, ".Rds"))

sim <- sim_huge_sf_p200[[42]]
which(sim$theta != 0, arr.ind = TRUE)

diagonals <- list(c(1,1), c(2,2), c(28,28), c(91,91), c(143,143), c(187,187))
non_zeros <- list(c(1,2), c(2,28), c(36,164), c(57,81), c(121,62), c(199, 132))
zeros <- list(c(1,80), c(3,10), c(17,50), c(45,100), c(109, 188), c(175, 180))

# Start of the trace plots
pdf("Figures/Supplementary/Huge_sf_data_42_traceplots.pdf", width = 7, 8.5)
cols <- palette.colors(6)
ltys <- rep(c(2,4,5,6,1), length.out = 6)
par(mfrow = c(3,1), mar = c(3.1, 3.1, 0.2, 0.2), mgp = c(2,1,0))
plot(huge_sf_dataset_42_p200[["42"]][1:1500, convert_array_inds_to_ltri_ind(diagonals[[1]])], type = "l",
     ylim = c(0, 8), lty = ltys[1], xlab = "Iteration", ylab = "Diagonal of the precision matrix")
for (i in 2:length(diagonals)) {
  lines(huge_sf_dataset_42_p200[["42"]][1:1500, convert_array_inds_to_ltri_ind(diagonals[[i]])],
        col = cols[i], lty = ltys[i])
}

plot(huge_sf_dataset_42_p200[["42"]][1:1500, convert_array_inds_to_ltri_ind(non_zeros[[1]])], type = "l",
     ylim = c(-0.5, 1.4), lty = ltys[1], xlab = "Iteration",
     ylab = "Non-zero off-diagonal of the precision matrix")
for (i in 2:length(non_zeros)) {
  lines(huge_sf_dataset_42_p200[["42"]][1:1500, convert_array_inds_to_ltri_ind(non_zeros[[i]])],
        col = cols[i], lty = ltys[i])
}

plot(huge_sf_dataset_42_p200[["42"]][1:1500, convert_array_inds_to_ltri_ind(zeros[[1]])], type = "l",
     ylim = c(-0.75, 0.75), lty = ltys[1], xlab = "Iteration",
     ylab = "Zero off-diagonal of the precision matrix")
for (i in 2:length(non_zeros)) {
  lines(huge_sf_dataset_42_p200[["42"]][1:1500,convert_array_inds_to_ltri_ind(zeros[[i]])],
        col = cols[i], lty = ltys[i])
}
dev.off()

# End of the trace plots
pdf("Figures/Supplementary/Huge_sf_data_42_traceplots_end.pdf", width = 7, 8.5)
cols <- palette.colors(6)
ltys <- rep(c(2,4,5,6,1), length.out = 6)
par(mfrow = c(3,1), mar = c(3.1, 3.1, 0.2, 0.2), mgp = c(2,1,0))
plot(1501:10000, huge_sf_dataset_42_p200[["42"]][1501:10000, convert_array_inds_to_ltri_ind(diagonals[[1]])],
     type = "l", ylim = c(0, 4), lty = ltys[1], xlab = "Iteration", ylab = "Diagonal of the precision matrix")
for (i in 2:length(diagonals)) {
  lines(1501:10000, huge_sf_dataset_42_p200[["42"]][1501:10000, convert_array_inds_to_ltri_ind(diagonals[[i]])],
        col = cols[i], lty = ltys[i])
}
plot(1501:10000, huge_sf_dataset_42_p200[["42"]][1501:10000, convert_array_inds_to_ltri_ind(non_zeros[[1]])],
     type = "l", ylim = c(-0.5, 1.4), lty = ltys[1], xlab = "Iteration",
     ylab = "Non-zero off-diagonal of the precision matrix")
for (i in 2:length(non_zeros)) {
  lines(1501:10000, huge_sf_dataset_42_p200[["42"]][1501:10000, convert_array_inds_to_ltri_ind(non_zeros[[i]])],
        col = cols[i], lty = ltys[i])
}
plot(1501:10000, huge_sf_dataset_42_p200[["42"]][1501:10000, convert_array_inds_to_ltri_ind(zeros[[1]])],
     type = "l", ylim = c(-0.75, 0.75), lty = ltys[1], xlab = "Iteration",
     ylab = "Zero off-diagonal of the precision matrix")
for (i in 2:length(non_zeros)) {
  lines(1501:10000, huge_sf_dataset_42_p200[["42"]][1501:10000, convert_array_inds_to_ltri_ind(zeros[[i]])],
        col = cols[i], lty = ltys[i])
}
dev.off()

# Remove datasets to save memory for further analysis.
rm(huge_sf_dataset_42_p200)
rm(sim_huge_sf_p200)
gc()

# Read MCMC samples from the MATLAB files.
huge_sf_5_datasets_p200 <- read_matlab_matrices(path, datanros, n_mcmc = 10000)

# Calculate posterior summaries.
huge_sf_post_summaries_p200 <- calculate_posterior_summaries(huge_sf_5_datasets_p200, datanros = datanros,
                                                             n_cores = 6)
# Combine scores over datasets.
huge_sf_combined_scores_p200 <- combine_over_datasets(huge_sf_post_summaries_p200, "rhat", datanros,
                                                      n_warmups, mcmc_lengths)

# Save only rhat, ess_bulk and ess_tail.
for (datanro in datanros) {
  for (n_warmup in n_warmups) {
    for (mcmc_length in mcmc_lengths) {
      huge_sf_post_summaries_p200[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(mcmc_length)]] <- huge_sf_post_summaries_p200[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(mcmc_length)]][,c("rhat", "ess_bulk", "ess_tail")]
    }
  }
}
for (n_warmup in n_warmups) {
  for (mcmc_length in mcmc_lengths) {
    huge_sf_combined_scores_p200[[paste0(n_warmup)]][[paste0(mcmc_length)]] <- huge_sf_combined_scores_p200[[paste0(n_warmup)]][[paste0(mcmc_length)]][,c("rhat", "ess_bulk", "ess_tail")]
  }
}

# Save summaries.
saveRDS(huge_sf_post_summaries_p200, file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/huge_sf_per_dataset_p200.Rds"))
saveRDS(huge_sf_combined_scores_p200, file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/huge_sf_combined_p200.Rds"))

# Read summaries if already done.
huge_sf_post_summaries_p200 <- readRDS(file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/huge_sf_per_dataset_p200.Rds"))
huge_sf_combined_scores_p200 <- readRDS(file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/huge_sf_combined_p200.Rds"))

# Plot per dataset boxplots for Rhats and ESSs.
par(mfrow = c(5,3))
for (datanro in datanros) {
  create_boxplots(huge_sf_post_summaries_p200[[paste0(datanro)]], p, "rhat", warm_ups = n_warmups,
                  mcmc_lengths = mcmc_lengths, ylim = c(1, 1.2), plot_one_row = FALSE)
}
par(mfrow = c(5,3))
for (datanro in datanros) {
  create_boxplots(huge_sf_post_summaries_p200[[paste0(datanro)]], p, "ess_bulk", warm_ups = n_warmups,
                  mcmc_lengths = mcmc_lengths, ylim = c(0, 7000), plot_one_row = FALSE)
}

# Create boxplots for combined results.
create_boxplots(huge_sf_combined_scores_p200, p, "rhat", n_warmups, mcmc_lengths, ylim = c(1, 1.2),
                filename = "Figures/Supplementary/huge_sf_p200_Rhats.png", image_size = c(11, 5),
                save_png = TRUE, res = 800)

create_boxplots(huge_sf_combined_scores_p200, p, "ess_bulk", n_warmups, mcmc_lengths, ylim = c(0, 7000),
                filename = "Figures/Supplementary/huge_sf_p200_ESSs.png", image_size = c(11, 5),
                save_png = TRUE, res = 800)

# Remove datasets to save memory for further analysis.
# rm(huge_sf_5_datasets_p200)
# gc()

### Not performed!!! MCMC diagnostics for hub network datasets generated using huge with p = 200.
# Initialize parameters.
# p <- 200
# structure <- "hubs"
# path_matrices <- paste0("Results_files/p", p, "/GHS_MCMC/",  structure,  "/")
# datanros <- c(9, 25, 33, 42, 47)
# n_warmups <- c(0, 500, 1000)
# mcmc_lengths <- seq(1000, 6000, 1000)
# 
# # Read MCMC samples from the MATLAB files.
# hub_5_datasets_p200 <- read_matlab_matrices(path, 9, n_mcmc = 10000)
# 
# # Calculate posterior summaries.
# hub_post_summaries_p200 <- calculate_posterior_summaries(hub_5_datasets_p200, datanros = datanros,
#                                                          n_cores = 16)
# # Combine scores over datasets.
# hub_combined_scores_p200 <- combine_over_datasets(hub_post_summaries_p200, "rhat", datanros,
#                                                   n_warmups, mcmc_lengths)
# 
# Save only rhat, ess_bulk and ess_tail.
# for (datanro in datanros) {
#   for (n_warmup in n_warmups) {
#     for (mcmc_length in mcmc_lengths) {
#       hub_post_summaries_p200[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(mcmc_length)]] <- hub_post_summaries_p200[[paste0(datanro)]][[paste0(n_warmup)]][[paste0(mcmc_length)]][,c("rhat", "ess_bulk", "ess_tail")]
#     }
#   }
# }
# for (n_warmup in n_warmups) {
#   for (mcmc_length in mcmc_lengths) {
#     hub_combined_scores_p200[[paste0(n_warmup)]][[paste0(mcmc_length)]] <- hub_combined_scores_p200[[paste0(n_warmup)]][[paste0(mcmc_length)]][,c("rhat", "ess_bulk", "ess_tail")]
#   }
# }
#
# # Save summaries.
# saveRDS(hub_post_summaries_p200, file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/hub_per_dataset.Rds"))
# saveRDS(hub_combined_scores_p200, file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/hub_combined_p200.Rds"))
# 
# # Read summaries if already done.
# hub_post_summaries_p200 <- readRDS(file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/hub_per_dataset_p200.Rds"))
# hub_combined_scores_p200 <- readRDS(file = paste0("Results_files/p", p, "/GHS_MCMC/MCMC_diagnostics/hub_combined_p200.Rds"))
# 
# # Plot per dataset boxplots for Rhats and ESSs.
# par(mfrow = c(5,3))
# for (datanro in datanros) {
#   create_boxplots(hub_post_summaries_p200[[paste0(datanro)]], p, "rhat", warm_ups = n_warmups,
#                   mcmc_lengths = mcmc_lengths, ylim = c(1, 1.2), plot_one_row = FALSE)
# }
# par(mfrow = c(5,3))
# for (datanro in datanros) {
#   create_boxplots(hub_post_summaries_p200[[paste0(datanro)]], p, "ess_bulk", warm_ups = n_warmups,
#                   mcmc_lengths = mcmc_lengths, ylim = c(0, 7000), plot_one_row = FALSE)
# }
# 
# # Create boxplots for combined results.
# create_boxplots(hub_combined_scores_p200, p, "rhat", n_warmups, mcmc_lengths, ylim = c(1, 1.2))
# 
# create_boxplots(hub_combined_scores_p200, p, "ess_bulk", n_warmups, mcmc_lengths, ylim = c(0, 7000))
# 
# Remove datasets to save memory for further analysis.
# rm(hub_5_datasets_p200)
# gc()
