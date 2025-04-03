# Load functions for plotting networks.
source("R_files/network_visualization.R")

# Set parameters.
n <- 120
p <- 100

# Set path (no needed to change).
path <- paste0("Data/n", n, "_p", p, "/")

# Load datasets.
sim_random <- readRDS(file = paste0(path, "bdgraph_random_n", n, "_p", p, ".Rds"))
sim_bdgraph_sf <- readRDS(file = paste0(path, "bdgraph_scale-free_n", n, "_p", p, ".Rds"))
sim_huge_sf <- readRDS(file = paste0(path, "huge_scale-free_n", n, "_p", p, ".Rds"))
sim_hubs <- readRDS(file = paste0(path, "huge_hubs_n", n, "_p", p, ".Rds"))

# Plot simulated network structures with p = 100 and save as eps-file.
setEPS()
postscript("Figures/Main_article/Fig_1_Simulated_structures_p100.eps", width = 20, height = 5)
par(mfrow = c(1,4))
set.seed(20250403)
plot_network(sim_random[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_bdgraph_sf[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_huge_sf[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_hubs[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
dev.off()

# Plot simulated network structures with p = 100 and save as pdf-file.
pdf("Figures/Supplementary/Simulated_structures_p100.pdf", width = 20, height = 5)
par(mfrow = c(1,4))
set.seed(20250403)
plot_network(sim_random[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_bdgraph_sf[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_huge_sf[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_hubs[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
dev.off()


# Set parameters.
n <- 120
p <- 200

# Set path (no needed to change).
path <- paste0("Data/n", n, "_p", p, "/")

# Load datasets.
sim_random <- readRDS(file = paste0(path, "bdgraph_random_n", n, "_p", p, ".Rds"))
sim_bdgraph_sf <- readRDS(file = paste0(path, "bdgraph_scale-free_n", n, "_p", p, ".Rds"))
sim_huge_sf <- readRDS(file = paste0(path, "huge_scale-free_n", n, "_p", p, ".Rds"))
sim_hubs <- readRDS(file = paste0(path, "huge_hubs_n", n, "_p", p, ".Rds"))

# Plot simulated network structures with p = 200 and save as pdf-file.
pdf("Figures/Supplementary/Simulated_structures_p200.pdf", width = 20, height = 5)
par(mfrow = c(1,4))
set.seed(20250403)
plot_network(sim_random[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_bdgraph_sf[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_huge_sf[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_hubs[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
dev.off()
