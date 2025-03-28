library(GHSGEM)

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

# Check if results file exists and run analysis if not.
if (!file.exists("Results_files/Analysis_of_inits/Hub_Thetas_100_inits.rds")) {
  # Do analysis using datasets with hub structure.
  cl <- makeCluster(detectCores(logical = FALSE))
  registerDoParallel(cl)
  (start <- Sys.time())
  Thetas_p100 <- foreach (data_nro = 1:n_datasets, .packages = "GHSGEM") %dopar% {
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
  saveRDS(Thetas_p100, file = "Results_files/Analysis_of_inits/Hub_Thetas_100_inits.rds")
}

# Load results if saved.
if (file.exists("Results_files/Analysis_of_inits/Hub_Thetas_100_inits.rds")) {
  Thetas_p100 <- readRDS("Results_files/Analysis_of_inits/Hub_Thetas_100_inits.rds")
}

# Check which elements (connections) are not exact zeros or ones.
indices <- c()
datanro <- c()
for (i in 1:n_datasets) {
  inds <- which(Thetas_p100[[i]] > 0 & Thetas_p100[[i]] < 1 & lower.tri(Thetas_p100[[i]]), arr.ind = TRUE)
  if (length(inds) > 0) {
    datanro <- c(datanro, i)
    indices <- c(indices, list(inds))
    print(paste0("Dataset nro ", i))
    print(inds)
    print(Thetas_p100[[i]][inds])
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

# Check if results file exists and run analysis if not.
if (!file.exists("Results_files/Analysis_of_inits/Huge_sf_Thetas_50_inits.rds")) {
  # Do analysis using datasets with scale-free structure (huge).
  cl <- makeCluster(detectCores(logical = FALSE))
  registerDoParallel(cl)
  (start <- Sys.time())
  Thetas_p200 <- foreach (data_nro = 1:n_datasets, .packages = "GHSGEM") %dopar% {
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
  saveRDS(Thetas_p200, file = "Results_files/Analysis_of_inits/Huge_sf_Thetas_50_inits.rds")
}

if (file.exists("Results_files/Analysis_of_inits/Huge_sf_Thetas_50_inits.rds")) {
  Thetas_p200 <- readRDS("Results_files/Analysis_of_inits/Huge_sf_Thetas_50_inits.rds")
}

# Check which elements (connections) are not exact zeros or ones.
indices <- c()
datanro <- c()
for (i in 1:n_datasets) {
  inds <- which(Thetas_p200[[i]] > 0 & Thetas_p200[[i]] < 1 & lower.tri(Thetas_p200[[i]]), arr.ind = TRUE)
  if (length(inds) > 0) {
    datanro <- c(datanro, i)
    indices <- c(indices, list(inds))
    print(paste0("Dataset nro ", i))
    print(inds)
    print(Thetas_p200[[i]][inds])
    print(paste0("Sum of those connections, which are not exact 0 or 1: ", length(inds)/2))
  }
}

# Check what is those connections estimates if default initial values are used.
for (i in 1:length(datanro)) {
  dn <- datanro[i]
  map <- GHS_MAP_estimation(sim_huge_sf[[dn]]$data, verbose = 0)
  print(map$Theta[indices[[i]]])
}
