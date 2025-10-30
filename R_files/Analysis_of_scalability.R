library(doParallel)

# Function which creates data and runs bots C++ and R version of algorithm, and calculates runtime.
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
plot(variables[1:5], mean_runtime_Cpp[1:5], type = "l", xlab = "Number of variables",
     ylab = "")
title(ylab = "Mean time per iteration (sec)", mgp = c(2.5, 1, 0))
grid()
text(100, mean_runtime_Cpp[5], "A", cex = 1.5)
# Plot full x-axis and both axes in logarithmic scale.
plot(variables, mean_runtime_Cpp, type = "l", log = "xy", xlab = "Number of variables",
     ylab = "")
title(ylab = "Mean time per iteration (sec)", mgp = c(2.5, 1, 0))
lines(c(100, 1200), mean_runtime_Cpp[c(1,9)], lty = 2)
grid()
text(100, mean_runtime_Cpp[9], "B", cex = 1.5)
options(scipen = 0)
dev.off()

