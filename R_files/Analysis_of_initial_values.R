library(GHSGEM)

generate_initial_values <- function(p, seed) {
  set.seed(seed)
  A <- matrix(rnorm(p^2, 0, 0.05), p, p)
  A <- (t(A) + A) / 2
  diag(A) <- 1
  return(A)
}

# Possible network structures.
structures <- c("random", "bdgraph_sf", "huge_sf", "hubs")

# Analysis for the datasets with 100 variables.
p <- 100

# Set path (no needed to change).
path <- paste0("Data/n", n, "_p", p, "/")

# Load datasets.
sim_hubs <- readRDS(file = paste0(path, "huge_hubs_n", n, "_p", p, ".Rds"))

# Set number of initial values.
n_inits <- 100
# Set number of datasets.
n_datasets <- length(sim_random)

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
if (file.exists("Results_files/Analysis_of_inits/Hubs_Thetas_100_inits.rda")) {
  Thetas_p100 <- readRDS("Results_files/Analysis_of_inits/Hubs_Thetas_100_inits.rda")
}

# Check which elements (connections) are not exact zeros or ones.
indices <- data.frame()
datanro <- c()
for (i in 1:n_data) {
  if (sum(which(Thetas_p100[[i]] > 0 & Thetas_p100[[i]] < 1)) > 0) {
    datanro <- c(datanro, i)
    indices <- rbind(indices, which(Thetas_p100[[i]] > 0 & Thetas_p100[[i]] < 1, arr.ind = TRUE)[1,])
    print(paste0("Dataset nro ", i))
    print(which(Thetas_p100[[i]] > 0 & Thetas_p100[[i]] < 1, arr.ind = TRUE))
    print(Thetas_p100[[i]][which(Thetas_p100[[i]] > 0 & Thetas_p100[[i]] < 1)])
    print(paste0("Sum of those connections, which are not exact 0 or 1: ", sum((Thetas_p100[[i]] > 0 & Thetas_p100[[i]] < 1))/2))
  }
}

# Check what is those connections estimates if default initial values are used.
for (i in 1:length(datanro)) {
  dn <- datanro[i]
  map <- GHS_MAP_estimation(sim_hubs[[dn]]$data, verbose = 0)
  print(map$Theta[indices[i,1], indices[i,2]])
}

########
# Analysis for the datasets with 200 variables.
p <- 200

# Load datasets.
sim_huge_sf <- readRDS(file = paste0(path, "huge_scale-free_n", n, "_p", p, ".Rds"))

# Set number of initial values.
n_inits <- 50
# Set number of datasets.
n_datasets <- length(sim_random)

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
indices <- data.frame()
datanro <- c()
for (i in 1:n_data) {
  if (sum(which(Thetas_p200[[i]] > 0 & Thetas_p200[[i]] < 1)) > 0) {
    datanro <- c(datanro, i)
    indices <- rbind(indices, which(Thetas_p200[[i]] > 0 & Thetas_p200[[i]] < 1, arr.ind = TRUE)[1,])
    print(paste0("Dataset nro ", i))
    print(which(Thetas_p100[[i]] > 0 & Thetas_p200[[i]] < 1, arr.ind = TRUE))
    print(Thetas_p100[[i]][which(Thetas_p200[[i]] > 0 & Thetas_p200[[i]] < 1)])
    print(paste0("Sum of those connections, which are not exact 0 or 1: ", sum((Thetas_p200[[i]] > 0 & Thetas_p200[[i]] < 1))/2))
  }
}

# Check what is those connections estimates if default initial values are used.
for (i in 1:length(datanro)) {
  dn <- datanro[i]
  map <- GHS_MAP_estimation(sim_huge_sf[[dn]]$data, verbose = 0)
  print(map$Theta[indices[i,1], indices[i,2]])
}



source("network_visualisation_and_scores.R")

generate_initial_values <- function(p, seed) {
  set.seed(seed)
  A <- matrix(runif(p^2)*2-1, ncol = p) 
  init <- round(cov2cor(t(A) %*% A), 10)
  if (matrixcalc::is.positive.definite(init)) {
    return(init)
  }
  else {
    cat(paste0("Seed ", seed, " does not produce positive definite matrix. Adding 1000 to the seed!"))
    generate_initial_values(p, seed + 1000)
  }
}

n <- 300
p <- 300

sim <- huge::huge.generator(n = n, d = p, graph = "scale-free")

initial_values1 <- array(dim = c(p, p, 50))
for (i in 1:50) {
  initial_values1[,,i] <- generate_initial_values(p, i + 20250210)
}

initial_values2 <- array(dim = c(p, p, 50))
for (i in 1:50) {
  initial_values2[,,i] <- generate_initial_values(p, i + 20260210)
}

Thetas1 <- Omegas1 <- Sigmas1 <- array(dim = c(p, p, 50))
for (i in 1:50) {
  map <- GHS_MAP_estimation(sim$data, use_Cpp = FALSE, initial_values = list(Omega = diag(1, p, p),
                                                                             Sigma = initial_values1[,,i]))
  Thetas1[,,i] <- map$Theta
  Omegas1[,,i] <- map$Omega_est
  Sigmas1[,,i] <- map$Sigma_est
}

sp1 <- matrix(50, p, p)
diag(sp1) <- 100

plot_matrix(sim$sigma, breaks = seq(-1, 1, 0.1))

map1 <- GHS_MAP_estimation(sim$data)

map2 <- GHS_MAP_estimation(sim$data, use_Cpp = FALSE, initial_values = list(Omega = matrix(100, p, p),
                                                                            Sigma = round(huge::huge.generator(n = n, d = p)$sigma), 12))

map3 <- GHS_MAP_estimation(sim$data, use_Cpp = FALSE, initial_values = list(Omega = round(sim$omega, 12),
                                                                            Sigma = round(sim$sigma, 12)))

Theta_mean <- apply(Thetas, 1:2, mean)

all(Theta_mean == Thetas[,,1])
all(Theta_mean == map1$Theta)
all(Theta_mean == map2$Theta)

sigma_F_norms <- c()
for (i in 1:50) {
  sigma_F_norms[i] <- print(norm(sim$sigma - Sigmas[,,i], type = "f"))
}
norm(sim$sigma - map2$Sigma_est, type = "f")
norm(sim$sigma, type = "f")

omega_F_norms <- c()
for (i in 1:50) {
  omega_F_norms[i] <- print(norm(sim$omega - Omegas[,,i], type = "f"))
}
norm(sim$omega - map2$Omega_est, type = "f")

sd(sigma_F_norms)
sd(omega_F_norms)

Theta_mean1 <- apply(Thetas1, 1:2, mean)

all(Theta_mean1 == Thetas1[,,1])
all(Theta_mean1 == map1$Theta)
all(Theta_mean1 == map2$Theta)



plot_matrix(initial_values2[,,5]*100)


generate_starting_points <- function(p, n_points, verbose = 0) {
  # Initialize Omega_saves matrix
  starting_points <- array(0, dim = c(p, p, n_points))
  
  for (i in 1:n_points) {
    # Initialize start_point as an identity matrix
    start_point <- diag(p)
    
    row <- 2
    
    for (row in 2:p) {
      d <- 0
      while (d != p) {
        row_seq <- row:p
        col_seq <- 1:(p - row + 1)
        
        #rand_noise <- -0.05 + runif(length(row_seq), 0, 1) * 2 * 0.05
        # Alternatively, use the following line for a different range:
        rand_noise <- -0.1 + runif(length(row_seq), 0, 1) * 2 * 0.1
        
        #lin_idcs <- as.vector(t(matrix(c(row_seq, col_seq), ncol = 2)))
        lin_idcs <- (col_seq - 1) * p + row_seq
        start_point[lin_idcs] <- rand_noise
        
        #lin_idcs <- as.vector(t(matrix(c(col_seq, row_seq), ncol = 2)))
        lin_idcs <- (row_seq - 1) * p + col_seq
        start_point[lin_idcs] <- rand_noise
        
        d <- sum(eigen(start_point)$values > 0)
      }
    }
    
    if (verbose > 0) {
      cat("Finished", i, "data set generation out of", n_points, "data sets\n")
    }
    
    # Store start_point in Omega_saves
    starting_points[,,i] <- start_point
  }
  return(starting_points)
}

spoints <- generate_starting_points(50, 1)[,,1]

par(mfrow = c(1,1))
plot_matrix(initial_values2[,,1])
plot_matrix(spoints)

hist(initial_values2[,,1][lower.tri(initial_values2[,,1])])
hist(spoints[lower.tri(spoints)])

plot_matrix(solve(spoints))
hist(solve(spoints)[lower.tri(solve(spoints))])

#######
library(GHSGEM)
library(doParallel)

generate_initial_values <- function(p, seed) {
  set.seed(seed)
  A <- matrix(runif(p^2)*2-1, ncol = p)
  init <- round(cov2cor(t(A) %*% A), 10)
  if (matrixcalc::is.positive.definite(init)) {
    return(init)
  }
  else {
    cat(paste0("Seed ", seed, " does not produce positive definite matrix. Adding 1000 to the seed!"))
    generate_initial_values(p, seed + 1000)
  }
}

generate_initial_values <- function(p, seed) {
  set.seed(seed)
  A <- matrix(rnorm(p^2, 0, 0.05), p, p)
  A <- (t(A) + A) / 2
  diag(A) <- 1
  return(A)
}

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

n_inits <- 100
n_data <- length(sim_random)
n_data <- 20

Sigma_inits <- Lambda_inits <- Delta_inits <- array(dim = c(p, p, n_inits))
for (i in 1:n_inits) {
  set.seed(20310218 + i)
  lambda <- matrix(rgamma(p^2, shape = 1, rate = 0.35), p, p)
  lambda <- (t(lambda) + lambda) / 2
  diag(lambda) <- 1
  Lambda_inits[,,i] <- lambda
  set.seed(20320218 + i)
  delta <- matrix(runif(p^2, 0, 100), p, p)
  delta <- (t(delta) + delta) / 2
  diag(delta) <- 1
  Delta_inits[,,i] <- delta
  Sigma_inits[,,i] <- generate_initial_values(p, 20300218 + i)
}

cl <- makeCluster(8)
registerDoParallel(cl)
start <- Sys.time()
Thetas <- foreach (data_nro = 1:n_data, .packages = "GHSGEM") %dopar% {
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
total_time <- Sys.time() - start
stopCluster(cl)
stopImplicitCluster()

for (i in 1:n_data) {
  print(paste0("Data nro ", i))
  print(which(Thetas[[i]] > 0 & Thetas[[i]] < 1))
  print(paste0("Nollasta ja ykkösestä poikkeavien määrä ", sum((Thetas[[i]] > 0 & Thetas[[i]] < 1))/2))
} 

#save(Thetas, file = "C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\Tuloksia\\Artikkeli\\p100\\Hubs_Thetas.rda")

load("C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\Tuloksia\\Artikkeli\\p100\\Hubs_Thetas.rda")

for (i in 1:20) {
  if (sum((Thetas[[i]] > 0 & Thetas[[i]] < 1)) > 0) {
    print(paste0("Data nro ", i))
    print(which(Thetas[[i]] > 0 & Thetas[[i]] < 1, arr.ind = TRUE))
    print(Thetas[[i]][which(Thetas[[i]] > 0 & Thetas[[i]] < 1, arr.ind = TRUE)])
    print(paste0("Nollasta ja ykkösestä poikkeavien määrä ", sum((Thetas[[i]] > 0 & Thetas[[i]] < 1))/2))
  }
}

Thetas[[15]][which(Thetas[[15]] > 0 & Thetas[[15]] < 1)]
Thetas[[20]][which(Thetas[[20]] > 0 & Thetas[[20]] < 1)]
Thetas[[21]][which(Thetas[[21]] > 0 & Thetas[[21]] < 1)]
Thetas[[29]][which(Thetas[[29]] > 0 & Thetas[[29]] < 1)]
Thetas[[32]][which(Thetas[[32]] > 0 & Thetas[[32]] < 1)]
Thetas[[41]][which(Thetas[[41]] > 0 & Thetas[[41]] < 1)]
Thetas[[46]][which(Thetas[[46]] > 0 & Thetas[[46]] < 1)]

sim_hubs[[15]]$sigma[75,73]
sim_hubs[[20]]$sigma[14,1]
sim_hubs[[21]]$sigma[5,1]
sim_hubs[[21]]$sigma[5,14]
sim_hubs[[29]]$sigma[35,25]
sim_hubs[[32]]$sigma[45,41]

inits <- generate_initial_values(50, 13448)
c2c_inits <- cov2cor(inits)
par(mfrow = c(2,2))
plot_matrix(inits, breaks = seq(-1, 1, 0.1))
hist(inits[lower.tri(inits)])
plot_matrix(c2c_inits, breaks = seq(-1, 1, 0.1))
hist(c2c_inits[lower.tri(c2c_inits)])

min(inits[lower.tri(inits)])
max(inits[lower.tri(inits)])

plot_matrix(sim_hubs[[10]]$sigmahat, breaks = seq(-1, 1, 0.1))

hist(lt_sigmahat[lt_sigma != 0])

lt_sigma <- sim_hubs[[10]]$sigma[lower.tri(sim_hubs[[10]]$sigma)]
lt_sigmahat <- sim_hubs[[10]]$sigmahat[lower.tri(sim_hubs[[10]]$sigmahat)]


for (i in 1:n_inits) {
  print(matrixcalc::is.positive.definite(Sigma_inits[,,i]))
}

for (i in 1:n_inits) {
  cat(min(Sigma_inits[,,i]), ",", max(Sigma_inits[,,i][lower.tri(Sigma_inits[,,i])]), "\n")
}

sim <- sim_hubs[[15]]

map <- GHS_MAP_estimation(sim$data, verbose = 1)

map$Theta[73, 75]

sim <- sim_hubs[[20]]

map <- GHS_MAP_estimation(sim$data, verbose = 1)

map$Theta[1, 14]

paste0("Start time: ", format(Sys.time(),'%H.%M.%S'))


?installed.packages
