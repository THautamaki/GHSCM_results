# Load functions for plotting networks.
source("R_files/network_visualization.R")

# Set parameter.
n <- 120
p <- 100

# Set path (no needed to change).
path <- paste0("Data/n", n, "_p", p, "/")

# Load datasets.
sim_random <- readRDS(file = paste0(path, "bdgraph_random_n", n, "_p", p, ".Rds"))
sim_bdgraph_sf <- readRDS(file = paste0(path, "bdgraph_scale-free_n", n, "_p", p, ".Rds"))
sim_huge_sf <- readRDS(file = paste0(path, "huge_scale-free_n", n, "_p", p, ".Rds"))
sim_hubs <- readRDS(file = paste0(path, "huge_hubs_n", n, "_p", p, ".Rds"))

# Plot simulated network structures and save as eps-file.
setEPS()
postscript("Figures/Main_article/Simulated_structures.eps", width = 20, height = 5)

par(mfrow = c(1,4))
plot_network(sim_random[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_bdgraph_sf[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_huge_sf[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_hubs[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)

dev.off()