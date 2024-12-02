library(huge)
library(BDgraph)

# Number of observations (sample size).
n <- 120
# Number of variables.
ps <- c(100, 200)
# Number of datasets.
n_datasets <- 50

for (p in ps) {
  # Folder where to save datasets. Change "path\\to\\folder" part where you want to save files.
  path <- paste0("C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\data\\Artikkeli\\p", p, "\\")
  #path <- paste0("path\\to\\folder\\p", p, "\\")
  
  # Initialise list where to save simulations.
  sim_random <- list()
  sim_bdgraph_sf <- list()
  sim_huge_sf <- list()
  sim_hubs <- list()
  
  for (i in 1:n_datasets) {
    # Set seed so same datasets can be generated again.
    seed <- 20241202 + i + p
    set.seed(seed)
    # Generate simulated dataset for random structure using bdgraph.sim function.
    sim <- bdgraph.sim(p = p, n = n, graph = "random", size = p/2)
    sim$omega <- sim$K
    sim$theta <- sim$G
    sim_random[[i]] <- sim
    
    # Set seed so same datasets can be generated again.
    seed <- 20251202 + i + p
    # Generate simulated dataset for random structure using bdgraph.sim function.
    sim <- bdgraph.sim(p = p, n = n, graph = "scale-free")
    sim$omega <- sim$K
    sim$theta <- sim$G
    sim_bdgraph_sf[[i]] <- sim
    
    # Set seed so same datasets can be generated again.
    seed <- 20261202 + i + p
    # Generate simulated dataset for random structure using bdgraph.sim function.
    sim <- huge.generator(n = n, d = p, graph = "scale-free")
    sim_huge_sf[[i]] <- sim
    
    # Set seed so same datasets can be generated again.
    seed <- 20271202 + i + p
    # Generate simulated dataset for random structure using bdgraph.sim function.
    sim <- huge.generator(n = n, d = p, graph = "hub")
    sim_hubs[[i]] <- sim
  }
  save(sim_random, file = paste0(path, "bdgraph_random_", p, ".Rda"))
  save(sim_bdgraph_sf, file = paste0(path, "bdgraph_scale-free_", p, ".Rda"))
  save(sim_huge_sf, file = paste0(path, "huge_scale-free_", p, ".Rda"))
  save(sim_hubs, file = paste0(path, "huge_hubs_", p, ".Rda"))
}

# Save necessary data also as csv-format for Matlab use.
for (p in ps) {
  # Path to folder where datasets are saved. Change "path\\to\\folder" part where you want to save and
  # load files.
  path <- paste0("C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\data\\Artikkeli\\p", p, "\\")
  #path <- paste0("path\\to\\folder\\p", p, "\\")
  
  # Load datasets.
  load(file = paste0(path, "bdgraph_random_", p, ".Rda"))
  load(file = paste0(path, "bdgraph_scale-free_", p, ".Rda"))
  load(file = paste0(path, "huge_scale-free_", p, ".Rda"))
  load(file = paste0(path, "huge_hubs_", p, ".Rda"))
  for (i in 1:n_datasets) {
    # Datasets with random structure from BDgraph.
    sim <- sim_random[[i]]
    write.csv(sim$data, paste0(path, "random\\random_data_nro_", i, ".csv"), row.names = FALSE)
    write.csv(sim$sigma, paste0(path, "random\\random_sigma_nro_", i, ".csv"), row.names = FALSE)
    write.csv(sim$omega, paste0(path, "random\\random_omega_nro_", i, ".csv"), row.names = FALSE)
    write.csv(sim$theta, paste0(path, "random\\random_theta_nro_", i, ".csv"), row.names = FALSE)
    
    # Datasets with scale-free structure from BDgraph.
    sim <- sim_bdgraph_sf[[i]]
    write.csv(sim$data, paste0(path, "bdgraph_scale-free\\bdgraph_sf_data_nro_", i, ".csv"), row.names = FALSE)
    write.csv(sim$sigma, paste0(path, "bdgraph_scale-free\\bdgraph_sf_sigma_nro_", i, ".csv"), row.names = FALSE)
    write.csv(sim$omega, paste0(path, "bdgraph_scale-free\\bdgraph_sf_omega_nro_", i, ".csv"), row.names = FALSE)
    write.csv(sim$theta, paste0(path, "bdgraph_scale-free\\bdgraph_sf_theta_nro_", i, ".csv"), row.names = FALSE)
    
    # Datasets with scale-free structure from huge.
    sim <- sim_huge_sf[[i]]
    write.csv(sim$data, paste0(path, "huge_scale-free\\huge_sf_data_nro_", i, ".csv"), row.names = FALSE)
    write.csv(sim$sigma, paste0(path, "huge_scale-free\\huge_sf_sigma_nro_", i, ".csv"), row.names = FALSE)
    write.csv(sim$omega, paste0(path, "huge_scale-free\\huge_sf_omega_nro_", i, ".csv"), row.names = FALSE)
    write.csv(as.matrix(sim$theta), paste0(path, "huge_scale-free\\huge_sf_theta_nro_", i, ".csv"), row.names = FALSE)
    
    # Datasets with hubs structure from huge.
    sim <- sim_hubs[[i]]
    write.csv(sim$data, paste0(path, "hubs\\hubs_data_nro_", i, ".csv"), row.names = FALSE)
    write.csv(sim$sigma, paste0(path, "hubs\\hubs_sigma_nro_", i, ".csv"), row.names = FALSE)
    write.csv(sim$omega, paste0(path, "hubs\\hubs_omega_nro_", i, ".csv"), row.names = FALSE)
    write.csv(as.matrix(sim$theta), paste0(path, "hubs\\hubs_theta_nro_", i, ".csv"), row.names = FALSE)
  }
}
