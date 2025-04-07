library(GHSGEM)
library(doParallel)

# Function which generates random positive-definite matrix.
generate_initial_values <- function(p, seed, sd = 0.1) {
  set.seed(seed)
  A <- matrix(rnorm(p^2, 0, sd), p, p)
  A <- (t(A) + A) / 2
  diag(A) <- 1 + max(abs(eigen(A, only.values = TRUE)$value))
  return(A)
}

# Possible network structures.
structures <- c("random", "bdgraph_sf", "huge_sf", "hubs")

# Analysis for the datasets with 100 variables.
n <- 120
p <- 100

# Set path (no needed to change).
path <- paste0("Data/n", n, "_p", p, "/")

# Load datasets.
sim_hubs <- readRDS(file = paste0(path, "huge_hubs_n", n, "_p", p, ".Rds"))

# Set number of initial values.
n_inits <- 100
# Set number of datasets.
n_datasets <- length(sim_hubs)

# Generate initial values.
Sigma_inits <- array(dim = c(p, p, n_inits))
for (i in 1:n_inits) {
  Sigma_inits[,,i] <- generate_initial_values(p, 20300218 + i)
}

hub_p100_resultsfile <- "Results_files/Analysis_of_inits/huge_hub_p100_Thetas_100_inits.rds"

# Check if results file exists and run analysis if not.
if (!file.exists(hub_p100_resultsfile)) {
  # Do analysis using datasets with hub structure.
  cl <- makeCluster(detectCores(logical = FALSE))
  registerDoParallel(cl)
  (start <- Sys.time())
  hub_Thetas_p100 <- foreach (data_nro = 1:n_datasets, .packages = "GHSGEM") %dopar% {
    theta <- matrix(0, p, p)
    for (i in 1:n_inits) {
      map <- GHS_MAP_estimation(sim_hubs[[data_nro]]$data, use_Cpp = FALSE,
                                initial_values = list(Omega = diag(1, p, p),
                                                      Sigma = Sigma_inits[,,i],
                                                      Lambda = matrix(1, p, p),
                                                      Delta = matrix(1, p, p)))
      theta <- theta + map$Theta
    }
    theta <- theta / n_inits
  }
  (total_time <- Sys.time() - start)
  stopCluster(cl)
  stopImplicitCluster()
  # Save results.
  saveRDS(hub_Thetas_p100, file = hub_p100_resultsfile)
}

# Load results if saved.
if (file.exists(hub_p100_resultsfile)) {
  hub_Thetas_p100 <- readRDS(hub_p100_resultsfile)
}

# Check which elements (connections) are not exact zeros or ones.
indices <- c()
datanro <- c()
for (i in 1:n_datasets) {
  inds <- which(hub_Thetas_p100[[i]] > 0 & hub_Thetas_p100[[i]] < 1 & lower.tri(hub_Thetas_p100[[i]]), arr.ind = TRUE)
  if (length(inds) > 0) {
    datanro <- c(datanro, i)
    indices <- c(indices, list(inds))
    print(paste0("Dataset nro ", i))
    print(inds)
    print(hub_Thetas_p100[[i]][inds])
    print(paste0("Sum of those connections, which are not exact 0 or 1: ", length(inds)/2))
  }
}

# Check what is those connections estimates if default initial values are used.
for (i in 1:length(datanro)) {
  dn <- datanro[i]
  map <- GHS_MAP_estimation(sim_hubs[[dn]]$data, verbose = 0)
  print(map$Theta[indices[[i]]])
}

########
# Analysis for the datasets with 200 variables.
n <- 120
p <- 200

# Set path (no needed to change).
path <- paste0("Data/n", n, "_p", p, "/")

# Load datasets.
sim_huge_sf <- readRDS(file = paste0(path, "huge_scale-free_n", n, "_p", p, ".Rds"))

# Set number of initial values.
n_inits <- 50
# Set number of datasets.
n_datasets <- length(sim_huge_sf)

# Generate initial values.
Sigma_inits <- array(dim = c(p, p, n_inits))
for (i in 1:n_inits) {
  Sigma_inits[,,i] <- generate_initial_values(p, 20300218 + i)
}

huge_sf_p200_resultsfile <- "Results_files/Analysis_of_inits/huge_sf_p200_Thetas_50_inits.rds"

# Check if results file exists and run analysis if not.
if (!file.exists(huge_sf_p200_resultsfile)) {
  # Do analysis using datasets with scale-free structure (huge).
  cl <- makeCluster(detectCores(logical = FALSE))
  registerDoParallel(cl)
  (start <- Sys.time())
  huge_sf_Thetas_p200 <- foreach (data_nro = 1:n_datasets, .packages = "GHSGEM") %dopar% {
    theta <- matrix(0, p, p)
    for (i in 1:n_inits) {
      map <- GHS_MAP_estimation(sim_huge_sf[[data_nro]]$data, use_Cpp = FALSE,
                                initial_values = list(Omega = diag(1, p, p),
                                                      Sigma = Sigma_inits[,,i],
                                                      Lambda = matrix(1, p, p),
                                                      Delta = matrix(1, p, p)))
      theta <- theta + map$Theta
    }
    theta <- theta / n_inits
  }
  (total_time <- Sys.time() - start)
  stopCluster(cl)
  stopImplicitCluster()
  # Save results.
  saveRDS(huge_sf_Thetas_p200, file = huge_sf_p200_resultsfile)
}

if (file.exists(huge_sf_p200_resultsfile)) {
  huge_sf_Thetas_p200 <- readRDS(huge_sf_p200_resultsfile)
}

# Check which elements (connections) are not exact zeros or ones.
indices <- c()
datanro <- c()
for (i in 1:n_datasets) {
  inds <- which(huge_sf_Thetas_p200[[i]] > 0 & huge_sf_Thetas_p200[[i]] < 1 & lower.tri(huge_sf_Thetas_p200[[i]]), arr.ind = TRUE)
  if (length(inds) > 0) {
    datanro <- c(datanro, i)
    indices <- c(indices, list(inds))
    print(paste0("Dataset nro ", i))
    print(inds)
    print(huge_sf_Thetas_p200[[i]][inds])
    print(paste0("Sum of those connections, which are not exact 0 or 1: ", length(inds)/2))
  }
}

# Check what is those connections estimates if default initial values are used.
for (i in 1:length(datanro)) {
  dn <- datanro[i]
  map <- GHS_MAP_estimation(sim_huge_sf[[dn]]$data, verbose = 0)
  print(map$Theta[indices[[i]]])
}


########
# Analysis for the datasets with 100 variables and BDgraph scale-free data.
n <- 120
p <- 100

# Set path (no needed to change).
path <- paste0("Data/n", n, "_p", p, "/")

# Load datasets.
sim_bdgraph_sf <- readRDS(file = paste0(path, "bdgraph_scale-free_n", n, "_p", p, ".Rds"))

# Set number of initial values.
n_inits <- 100
# Set number of datasets.
n_datasets <- length(sim_bdgraph_sf)

# Generate initial values.
Sigma_inits <- array(dim = c(p, p, n_inits))
for (i in 1:n_inits) {
  Sigma_inits[,,i] <- generate_initial_values(p, 20300218 + i)
}

bdgraph_sf_p100_resultsfile <- "Results_files/Analysis_of_inits/BDgraph_sf_p100_Thetas_100_inits.rds"

# Check if results file exists and run analysis if not.
if (!file.exists(bdgraph_sf_p100_resultsfile)) {
  # Do analysis using datasets with scale-free structure (huge).
  cl <- makeCluster(detectCores(logical = FALSE)-1)
  registerDoParallel(cl)
  print(start <- Sys.time())
  BDgraph_sf_Thetas_p100 <- foreach (data_nro = 1:n_datasets, .packages = "GHSGEM") %dopar% {
    theta <- matrix(0, p, p)
    for (i in 1:n_inits) {
      map <- GHS_MAP_estimation(sim_bdgraph_sf[[data_nro]]$data, use_Cpp = FALSE,
                                initial_values = list(Omega = diag(1, p, p),
                                                      Sigma = Sigma_inits[,,i],
                                                      Lambda = matrix(1, p, p),
                                                      Delta = matrix(1, p, p)))
      theta <- theta + map$Theta
    }
    theta <- theta / n_inits
  }
  print(total_time <- Sys.time() - start)
  stopCluster(cl)
  stopImplicitCluster()
  
  # Save results.
  saveRDS(BDgraph_sf_Thetas_p100, file = bdgraph_sf_p100_resultsfile)
}

if (file.exists(bdgraph_sf_p100_resultsfile)) {
  BDgraph_sf_Thetas_p100 <- readRDS(bdgraph_sf_p100_resultsfile)
}

# Check which elements (connections) are not exact zeros or ones.
indices <- datanro <- n_conns <- c()
for (i in 1:n_datasets) {
  inds <- which(BDgraph_sf_Thetas_p100[[i]] > 0 & BDgraph_sf_Thetas_p100[[i]] < 1 & lower.tri(BDgraph_sf_Thetas_p100[[i]]), arr.ind = TRUE)
  if (length(inds) > 0) {
    datanro <- c(datanro, i)
    indices <- c(indices, list(inds))
    n_conns <- c(n_conns, length(inds) / 2)
    print(paste0("Dataset nro ", i))
    print(paste0("Sum of those connections, which are not exact 0 or 1: ", length(inds)/2))
  }
}

# Check number of datasets which have differences in estimated connections.
length(datanro)

# Print ordered number of connections and calculate how many are under 10, between 10 and 30
# and over 30.
n_conns[order(n_conns)]
length(which(n_conns < 10))
length(which(n_conns >= 10 & n_conns <= 30))
length(which(n_conns > 30))

# Calculate the estimates using the default initial value.
maps <- c()
for (i in 1:n_datasets) {
  map <- GHS_MAP_estimation(sim_bdgraph_sf[[i]]$data, verbose = 0)
  maps <- c(maps, list(map))
}

# Check those which should be same as using identity matrix.
all_zero_or_one <- (1:50)[-datanro]
for (dn in all_zero_or_one) {
  cat("Data nro", dn, "\n")
  theta_est <- BDgraph_sf_Thetas_p100[[dn]]
  cm <- conf_matrix(sim_bdgraph_sf[[dn]]$theta, theta_est)
  cm2 <- conf_matrix(sim_bdgraph_sf[[dn]]$theta, maps[[dn]]$Theta_est)
  print(conf_matrix(maps[[dn]]$Theta_est, theta_est))
}

# Calculate scores over all datasets. Set all non-zero elements of the adjacency matrix to 1.
scores_100_inits <- scores_default <- data.frame()
for (dn in 1:n_datasets) {
  theta_est <- BDgraph_sf_Thetas_p100[[dn]]
  theta_est[theta_est != 0] <- 1
  cm <- conf_matrix(sim_bdgraph_sf[[dn]]$theta, theta_est)
  cm2 <- conf_matrix(sim_bdgraph_sf[[dn]]$theta, maps[[dn]]$Theta_est)
  scores_100_inits <- rbind(scores_100_inits, calculate_scores(cm))
  scores_default <- rbind(scores_default, calculate_scores(cm2))
}

# Set scores.
scores <- c("MCC", "F1", "TPR", "FPR", "FDR")

# Print mean of the scores.
round(colMeans(scores_default[, scores]), 3)
round(colMeans(scores_100_inits[, scores]), 3)
# Print sd of the scores.
round(apply(scores_default[, scores], 2, sd), 3)
round(apply(scores_100_inits[, scores], 2, sd), 3)

# Print scores for those datasets which differences are 30 connections or over.
round(scores_default[datanro[which(n_conns >= 30)], scores], 4)
round(scores_100_inits[datanro[which(n_conns >= 30)], scores], 4)

