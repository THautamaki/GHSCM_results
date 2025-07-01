# Load functions for plotting networks.
source("R_files/network_visualization.R")

# Load datasets.
# Set parameters.
n <- 120
p <- 100

# Set path (no needed to change).
path <- paste0("Data/n", n, "_p", p, "/")

# Load datasets.
sim_random_p100 <- readRDS(file = paste0(path, "bdgraph_random_n", n, "_p", p, ".Rds"))
sim_bdgraph_sf_p100 <- readRDS(file = paste0(path, "bdgraph_scale-free_n", n, "_p", p, ".Rds"))
sim_huge_sf_p100 <- readRDS(file = paste0(path, "huge_scale-free_n", n, "_p", p, ".Rds"))
sim_hubs_p100 <- readRDS(file = paste0(path, "huge_hubs_n", n, "_p", p, ".Rds"))

# Set parameters.
p <- 200

# Set path (no needed to change).
path <- paste0("Data/n", n, "_p", p, "/")

# Load datasets.
sim_random_p200 <- readRDS(file = paste0(path, "bdgraph_random_n", n, "_p", p, ".Rds"))
sim_bdgraph_sf_p200 <- readRDS(file = paste0(path, "bdgraph_scale-free_n", n, "_p", p, ".Rds"))
sim_huge_sf_p200 <- readRDS(file = paste0(path, "huge_scale-free_n", n, "_p", p, ".Rds"))
sim_hubs_p200 <- readRDS(file = paste0(path, "huge_hubs_n", n, "_p", p, ".Rds"))

# Plot simulated network structures with p = 100 and save as eps-file.
setEPS()
postscript("Figures/Main_article/Simulated_structures_p100.eps", width = 20, height = 5)
par(mfrow = c(1,4))
set.seed(20250403)
plot_network(sim_random_p100[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_bdgraph_sf_p100[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_huge_sf_p100[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_hubs_p100[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
dev.off()

# Plot simulated network structures with p = 100 and save as pdf-file.
pdf("Figures/Supplementary/Simulated_structures_p100.pdf", width = 20, height = 5)
par(mfrow = c(1,4))
set.seed(20250403)
plot_network(sim_random_p100[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_bdgraph_sf_p100[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_huge_sf_p100[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_hubs_p100[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
dev.off()


# Plot simulated network structures with p = 200 and save as eps-file.
setEPS()
postscript("Figures/Main_article/Simulated_structures_p200.eps", width = 20, height = 5)
par(mfrow = c(1,4))
set.seed(20250403)
plot_network(sim_random_p200[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
title(main = "A", adj = 0)
plot_network(sim_bdgraph_sf_p200[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_huge_sf_p200[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_hubs_p200[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)

dev.off()

#######


# Plot simulated network structures with both sizes in one figure and save as pdf-file.
pdf("Figures/Supplementary/Simulated_structures.pdf", width = 20, height = 10)
par(mfrow = c(2,4))
set.seed(20250403)
plot_network(sim_random_p100[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
title(main = "A", adj = 0.01, line = -2, cex.main = 3)
plot_network(sim_bdgraph_sf_p100[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_huge_sf_p100[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_hubs_p100[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
set.seed(20250403)
plot_network(sim_random_p200[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
title(main = "B", adj = 0.01, line = -2, cex.main = 3)
plot_network(sim_bdgraph_sf_p200[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_huge_sf_p200[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_hubs_p200[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
dev.off()

# Plot simulated network structures with both sizes in one figure and save as eps-file.
setEPS()
postscript("Figures/Main_article/Simulated_structures.eps", width = 20, height = 10)
par(mfrow = c(2,4))
set.seed(20250403)
plot_network(sim_random_p100[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
title(main = "A", adj = 0.01, line = -2, cex.main = 3)
plot_network(sim_bdgraph_sf_p100[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_huge_sf_p100[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_hubs_p100[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
set.seed(20250403)
plot_network(sim_random_p200[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
title(main = "B", adj = 0.01, line = -2, cex.main = 3)
plot_network(sim_bdgraph_sf_p200[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_huge_sf_p200[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_hubs_p200[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
dev.off()

