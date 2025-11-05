generate_datasets <- function(sample_sizes = 120, variable_numbers = c(100, 200), n_datasets = 50,
                              save_datasets = TRUE, return = FALSE) {
  # Initialize list for return.
  simulations <- list()
  for (n in sample_sizes) {
    for (p in variable_numbers) {
      cat(paste0("Creating datasets with dimensions n = ", n, " and p = ", p, ".\n"))
      # Folders where to save datasets. The default is to save these /Data/nxxx_pxxx/structure folder
      # into the working directory. The function creates these, if they not exists.
      if(!file.exists("Data")) {
        dir.create(file.path(getwd(), "Data"))
      }
      path <- paste0("Data/n", n, "_p", p)
      if (!file.exists(path)) {
        dir.create(file.path(getwd(), path))
      }
      if (!file.exists(paste0(path, "/random"))) {
        dir.create(file.path(getwd(), path, "random"))
      }
      if (!file.exists(paste0(path, "/bdgraph_sf"))) {
        dir.create(file.path(getwd(), path, "bdgraph_sf"))
      }
      if (!file.exists(paste0(path, "/huge_sf"))) {
        dir.create(file.path(getwd(), path, "huge_sf"))
      }
      if (!file.exists(paste0(path, "/hubs"))) {
        dir.create(file.path(getwd(), path, "hubs"))
      }
      # Initialise list where to save simulations.
      simulations[[paste0("n", n, "_p", p)]][["random"]]
      simulations[[paste0("n", n, "_p", p)]][["bdgraph_sf"]]
      simulations[[paste0("n", n, "_p", p)]][["huge_sf"]]
      simulations[[paste0("n", n, "_p", p)]][["hubs"]]
      
      for (i in 1:n_datasets) {
        cat(paste0("Creating dataset nro ", i, ".\n"))
        # Set seed so the same datasets can be generated again.
        seed <- 20241202 + i + p
        set.seed(seed)
        # Generate simulated dataset with random structure using bdgraph.sim function.
        sim1 <- BDgraph::bdgraph.sim(p = p, n = n, graph = "random", size = floor(p/2))
        sim1$omega <- sim1$K
        sim1$theta <- sim1$G
        simulations[[paste0("n", n, "_p", p)]][["random"]][[i]] <- sim1
        
        # Generate simulated dataset with random structure using bdgraph.sim function.
        sim2 <- BDgraph::bdgraph.sim(p = p, n = n, graph = "scale-free")
        sim2$omega <- sim2$K
        sim2$theta <- sim2$G
        simulations[[paste0("n", n, "_p", p)]][["bdgraph_sf"]][[i]] <- sim2
        
        # Generate simulated dataset with random structure using bdgraph.sim function.
        sim3 <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = FALSE)
        simulations[[paste0("n", n, "_p", p)]][["huge_sf"]][[i]] <- sim3
        
        # Generate simulated dataset with random structure using bdgraph.sim function.
        sim4 <- huge::huge.generator(n = n, d = p, graph = "hub", verbose = FALSE)
        simulations[[paste0("n", n, "_p", p)]][["hubs"]][[i]] <- sim4
      }
      
      # If save_datasets = TRUE, save they as Rda files and csv files for MATLAB use.
      if (save_datasets) {
        cat("Saving datasets...")
        sim_random <- simulations[[paste0("n", n, "_p", p)]][["random"]]
        saveRDS(sim_random, file = paste0(path, "/bdgraph_random_n", n, "_p", p, ".Rds"))
        sim_bdgraph_sf <- simulations[[paste0("n", n, "_p", p)]][["bdgraph_sf"]]
        saveRDS(sim_bdgraph_sf, file = paste0(path, "/bdgraph_scale-free_n", n, "_p", p, ".Rds"))
        sim_huge_sf <- simulations[[paste0("n", n, "_p", p)]][["huge_sf"]]
        saveRDS(sim_huge_sf, file = paste0(path, "/huge_scale-free_n", n, "_p", p, ".Rds"))
        sim_hubs <- simulations[[paste0("n", n, "_p", p)]][["hubs"]]
        saveRDS(sim_hubs, file = paste0(path, "/huge_hubs_n", n, "_p", p, ".Rds"))
        
        for (i in 1:n_datasets) {
          sim <- simulations[[paste0("n", n, "_p", p)]][["random"]][[i]]
          write.csv(sim$data, paste0(path, structure, "/", structure, "_data_nro_", i, ".csv"), row.names = FALSE)
          write.csv(sim$sigma, paste0(path, structure, "/", structure, "_sigma_nro_", i, ".csv"), row.names = FALSE)
          write.csv(sim$omega, paste0(path, structure, "/", structure, "_omega_nro_", i, ".csv"), row.names = FALSE)
          write.csv(sim$theta, paste0(path, structure, "/", structure, "_theta_nro_", i, ".csv"), row.names = FALSE)
          
          # Datasets with scale-free structure from BDgraph.
          sim <-simulations[[paste0("n", n, "_p", p)]][["bdgraph_sf"]][[i]]
          write.csv(sim$data, paste0(path, structure, "/", structure, "_data_nro_", i, ".csv"), row.names = FALSE)
          write.csv(sim$sigma, paste0(path, structure, "/", structure, "_sigma_nro_", i, ".csv"), row.names = FALSE)
          write.csv(sim$omega, paste0(path, structure, "/", structure, "_omega_nro_", i, ".csv"), row.names = FALSE)
          write.csv(sim$theta, paste0(path, structure, "/", structure, "_theta_nro_", i, ".csv"), row.names = FALSE)
          
          # Datasets with scale-free structure from huge.
          sim <- simulations[[paste0("n", n, "_p", p)]][["huge_sf"]][[i]]
          write.csv(sim$data, paste0(path, structure, "/", structure, "_data_nro_", i, ".csv"), row.names = FALSE)
          write.csv(sim$sigma, paste0(path, structure, "/", structure, "_sigma_nro_", i, ".csv"), row.names = FALSE)
          write.csv(sim$omega, paste0(path, structure, "/", structure, "_omega_nro_", i, ".csv"), row.names = FALSE)
          write.csv(as.matrix(sim$theta), paste0(path, structure, "/", structure, "_theta_nro_", i, ".csv"), row.names = FALSE)
          
          # Datasets with hubs structure from huge.
          sim <- simulations[[paste0("n", n, "_p", p)]][["hubs"]][[i]]
          write.csv(sim$data, paste0(path, structure, "/", structure, "_data_nro_", i, ".csv"), row.names = FALSE)
          write.csv(sim$sigma, paste0(path, structure, "/", structure, "_sigma_nro_", i, ".csv"), row.names = FALSE)
          write.csv(sim$omega, paste0(path, structure, "/", structure, "_omega_nro_", i, ".csv"), row.names = FALSE)
          write.csv(as.matrix(sim$theta), paste0(path, structure, "/", structure, "_theta_nro_", i, ".csv"), row.names = FALSE)
        }
        cat("done!\n")
      }
    }
  }
  if (return) {
    return(simulations)
  }
}

