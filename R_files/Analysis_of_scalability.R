library(doParallel)

# Function which creates data and runs both C++ and R version of algorithm, and calculates wall time.
run_one_p_parallel <- function(p, n = 0, n_repeats = 10, seed = NULL, n_iters = 20,
                               structure = "scale-free", n_cores = 0, verbose = 0) {
  if (n_cores == 0) n_cores <- min(n_repeats, detectCores() - 1)
  cat("p = ", p)
  n <- 0.6 * p
  set.seed(seed)
  sim <- huge::huge.generator(n = n, d = p, graph = structure, verbose = FALSE)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  runtimes_Cpp <- foreach(i = 1:n_repeats, .packages = "GHSCM") %dopar% {
    start <- Sys.time()
    map <- GHSCM::GHS_MAP_estimation(sim$data, max_iter = n_iters, verbose = verbose)
    end <- Sys.time()
    as.numeric(difftime(end, start, units = "secs"))
  }
  stopCluster(cl)
  stopImplicitCluster()
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  runtimes_R <- foreach(i = 1:n_repeats, .packages = "GHSCM") %dopar% {
    start <- Sys.time()
    map <- GHSCM::GHS_MAP_estimation(sim$data, max_iter = n_iters, verbose = verbose,
                                     use_Cpp = FALSE)
    end <- Sys.time()
    as.numeric(difftime(end, start, units = "secs"))
  }
  stopCluster(cl)
  stopImplicitCluster()
  return(list(runtimes_Cpp = runtimes_Cpp, runtimes_R = runtimes_R))
}

# Function which creates data, runs algorithm to full convergence and calculates wall time.
run_one_p_parallel_full <- function(p, np_ratios = 0, n_repeats = 10, seed = NULL, n_iters = 20,
                                    structure = "scale-free", n_cores = 0, verbose = 0) {
  if (n_cores == 0) n_cores <- min(n_repeats, detectCores() - 1)
  cat("p =", p, "\n")
  ns <- round(np_ratios * p)
  all_runtimes <- data.frame()
  set.seed(seed)
  for (n in ns) {
    cat("n =", n, "\n")
    sims <- list()
    for (i in 1:n_repeats) {
      sims[[i]] <- huge::huge.generator(n = n, d = p, graph = structure, verbose = FALSE)
    }
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    runtimes <- foreach(i = 1:n_repeats, .packages = "GHSCM", .combine = "rbind") %dopar% {
      start <- Sys.time()
      map <- GHSCM::GHS_MAP_estimation(sims[[i]]$data, max_iter = n_iters, verbose = verbose)
      end <- Sys.time()
      time <- as.numeric(difftime(end, start, units = "secs"))
      cbind(p = p, n = n, time = time, n_iters = map$iters)
    }
    stopCluster(cl)
    stopImplicitCluster()
    all_runtimes <- rbind(all_runtimes, runtimes)
  }
  return(all_runtimes)
}

#####
### Per iteration run time analysis with different number of variables.
# Define seeds for data creation, number of variables, number of iterations, and number of replications.
seeds <- c(18102025, 19102025, 20102025, 21102025, 28102025, 22102025, 23102025, 24102025, 25102025)
variables <- c(100, 200, 300, 400, 500, 600, 800, 1000, 1200)
n_iters <- 10
n_repeats <- 10

# Run analysis if results files do not exist.
if (!file.exists("Results_files/Runtimes/Runtimes_Cpp_p100-1200.rds")) {
  all_runtimes_Cpp <- all_runtimes_R <- data.frame()
  for (i in 1:length(variables)) {
    p <- variables[i]
    runtimes <- run_one_p_parallel(p, n_repeats = n_repeats, seed = seeds[i], n_iters = n_iters, verbose = 0)
    all_runtimes_Cpp <- rbind(all_runtimes_Cpp, runtimes$runtimes_Cpp)
    all_runtimes_R <- rbind(all_runtimes_R, runtimes$runtimes_R)
  }
  saveRDS(all_runtimes_Cpp, file = "Results_files/Runtimes/Runtimes_Cpp_p100-1200.rds")
  saveRDS(all_runtimes_R, file = "Results_files/Runtimes/Runtimes_R_p100-1200.rds")
}

# Read results files if they exist.
if (file.exists("Results_files/Runtimes/Runtimes_Cpp_p100-1200.rds")) {
  all_runtimes_Cpp <- readRDS("Results_files/Runtimes/Runtimes_Cpp_p100-1200.rds")
  all_runtimes_R <- readRDS("Results_files/Runtimes/Runtimes_R_p100-1200.rds")
}

# Calculate mean runtimes per iteration.
mean_runtime_Cpp <- rowMeans(all_runtimes_Cpp) / n_iters
mean_runtime_R <- rowMeans(all_runtimes_R) / n_iters

# Plot runtimes and save as EPS image.
setEPS()
postscript("Figures/Main_article/runtime_scaling.eps", width = 10, height = 4.6)
options(scipen = 999)
par(mfrow = c(1,2))
par(mgp = c(2, 1, 0))
par(mar = c(3.1, 4.1, 0.4, 0.2))
# Plot zoomed x-axis (100-500 variables) in normal scale.
plot(variables[1:5], mean_runtime_Cpp[1:5], type = "l", xlab = "Number of variables p",
     ylab = "")
title(ylab = "Per iteration wall time (sec)", mgp = c(2.5, 1, 0))
grid()
text(100, mean_runtime_Cpp[5], "A", cex = 1.5)
# Plot full x-axis and both axes in logarithmic scale.
plot(variables, mean_runtime_Cpp, type = "l", log = "xy", xlab = "Number of variables p",
     ylab = "")
title(ylab = "Per iteration wall time (sec)", mgp = c(2.5, 1, 0))
lines(c(100, 1200), mean_runtime_Cpp[c(1,9)], lty = 2)
grid()
text(100, mean_runtime_Cpp[9], "B", cex = 1.5)
options(scipen = 0)
dev.off()

#####
### Full run time analysis with different number of variables.
# Define seeds for data creation, number of variables, number of iterations, and number of replications.
seeds <- c(18102025, 19102025, 20102025, 21102025, 28102025, 22102025, 23102025, 24102025, 25102025)
variables <- c(100, 200, 300, 400, 500, 600, 800, 1000, 1200)
n_iters <- 1000
n_repeats <- 10
np_ratios <- c(0.6, 0.9, 1.2, 1.5)

# Run analysis if results files do not exist.
if (!file.exists("Results_files/Runtimes/Full_runtimes_p100-1200.rds")) {
  all_runtimes <- data.frame()
  for (i in 1:length(variables)) {
    p <- variables[i]
    runtimes <- run_one_p_parallel_full(p, np_ratios = np_ratios, n_repeats = n_repeats,
                                        seed = seeds[i], n_iters = n_iters, verbose = 0)
    all_runtimes <- rbind(all_runtimes, runtimes)
  }
  saveRDS(all_runtimes, file = "Results_files/Runtimes/Full_runtimes_p100-1200.rds")
}

# Read results file if it exists.
if (file.exists("Results_files/Runtimes/Full_runtimes_p100-1200.rds")) {
  all_runtimes <- readRDS("Results_files/Runtimes/Full_runtimes_p100-1200.rds")
}

# Aggregate results by n and p.
runtimes_agg <- aggregate(all_runtimes$time, by = list(all_runtimes$n, all_runtimes$p), FUN = mean)
runtimes_agg$np_ratio <- runtimes_agg$Group.1 / runtimes_agg$Group.2
iterations_agg <- aggregate(all_runtimes$n_iters, by = list(all_runtimes$n, all_runtimes$p), FUN = mean)
iterations_agg$np_ratio <- iterations_agg$Group.1 / runtimes_agg$Group.2

#####
# Create body of LaTeX table (Table 6 in the manuscript).
run <- TRUE
if (run) {
  cat("One iteration ")
  for (i in 1:length(mean_runtime_Cpp)) {
    cat(" & ", round(mean_runtime_Cpp[i], as.numeric(select_rounding(mean_runtime_Cpp[i])[1])))
  }
  cat(" \\\\ \n")
  for (j in 1:4) {
    cat("$n/p$ ratio ", np_ratios[j])
    times <- runtimes_agg[runtimes_agg$np_ratio == np_ratios[j],]
    for (i in 1:length(times$x)) {
      cat(" & ", round(times$x[i], as.numeric(select_rounding(times$x[i])[1])))
    }
    cat(" \\\\ \n")
  }
}

# Calculate empirical computational complexity.
(log(mean_runtime_Cpp[9]) - log(mean_runtime_Cpp[1])) / (log(1200) - log(100))
(log(runtimes_agg[runtimes_agg$np_ratio == np_ratios[1],][9, "x"]) - log(runtimes_agg[runtimes_agg$np_ratio == np_ratios[1],][1, "x"])) / (log(1200) - log(100))
(log(runtimes_agg[runtimes_agg$np_ratio == np_ratios[2],][9, "x"]) - log(runtimes_agg[runtimes_agg$np_ratio == np_ratios[2],][1, "x"])) / (log(1200) - log(100))
(log(runtimes_agg[runtimes_agg$np_ratio == np_ratios[3],][9, "x"]) - log(runtimes_agg[runtimes_agg$np_ratio == np_ratios[3],][1, "x"])) / (log(1200) - log(100))
(log(runtimes_agg[runtimes_agg$np_ratio == np_ratios[4],][9, "x"]) - log(runtimes_agg[runtimes_agg$np_ratio == np_ratios[4],][1, "x"])) / (log(1200) - log(100))


#####
# Create and save run time figure (not used in the manuscript).
# setEPS()
# postscript("Figures/Main_article/runtime_scaling.eps", width = 10.5, height = 3.6)
# options(scipen = 999)
# par(mfrow = c(1,3))
# par(mgp = c(2, 1, 0))
# par(mar = c(3.1, 3.1, 0.7, 0.2))
# # Plot zoomed x-axis (100-500 variables) in normal scale.
# plot(variables[1:5], mean_runtime_Cpp[1:5], type = "l", xlab = "Number of variables p",
#      ylab = "")
# title(ylab = "Wall time for one iteration (sec)")
# grid()
# text(100, mean_runtime_Cpp[5], "A", cex = 1.5)
# # Plot full x-axis and both axes in logarithmic scale.
# plot(variables, mean_runtime_Cpp, type = "l", log = "xy", xlab = "Number of variables p",
#      ylab = "")
# title(ylab = "Wall time for one iteration (sec)")
# lines(c(100, 1200), mean_runtime_Cpp[c(1,9)], lty = 2)
# grid()
# text(100, mean_runtime_Cpp[9], "B", cex = 1.5)
# 
# ltys <- c(1,4,5,6)
# colors <- palette.colors(4)
# options(scipen = 999)
# plot(runtimes_agg$Group.2[runtimes_agg$np_ratio == 0.6], runtimes_agg$x[runtimes_agg$np_ratio == 0.6]/60,
#      type = "l", log = "xy", lty = ltys[1], col = colors[1], xlab = "Number of variables p",
#      ylab = "Wall time to reach convergence (min)", ylim = c(min(runtimes_agg$x)/60, max(runtimes_agg$x)/60))
# lines(runtimes_agg$Group.2[runtimes_agg$np_ratio == 0.9], runtimes_agg$x[runtimes_agg$np_ratio == 0.9]/60,
#       lty = ltys[2], col = colors[2])
# lines(runtimes_agg$Group.2[runtimes_agg$np_ratio == 1.2], runtimes_agg$x[runtimes_agg$np_ratio == 1.2]/60,
#       lty = ltys[3], col = colors[3])
# lines(runtimes_agg$Group.2[runtimes_agg$np_ratio == 1.5], runtimes_agg$x[runtimes_agg$np_ratio == 1.5]/60,
#       lty = ltys[4], col = colors[4])
# grid()
# legend("bottomright", legend = c("n/p ratio = 0.6", "n/p ratio = 0.9", "n/p ratio = 1.2", "n/p ratio = 1.5"),
#        lty = ltys, col = colors, bg = "white")
# text(100, max(runtimes_agg$x)/60, "C", cex = 1.5)
# options(scipen = 0)
# dev.off()
