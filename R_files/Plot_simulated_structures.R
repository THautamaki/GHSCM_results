library(igraph)

plot_network <- function(a_mat, title = "", vertex_size = 3, margins = c(5.1, 4.1, 4.1, 2.1),
                         delete_isolates = FALSE, delete_nodes_degree = 0,
                         node_label_dist = 0.8, node_labels = NULL) {
  par(mar = margins)
  network <- graph_from_adjacency_matrix(a_mat, mode = "undirected", diag = F)
  if (delete_isolates) {
    network <- delete_vertices(network, V(network)[degree(network) == 0])
  }
  if (delete_nodes_degree > 0) {
    network <- delete_vertices(network, V(network)[degree(network) <= delete_nodes_degree])
  }
  coords <- layout_with_fr(network)
  if (is.null(node_labels)) {
    plot(network, layout = coords, vertex.label.color = "black", edge.width = 3,
         vertex.label.dist = node_label_dist, vertex.size = vertex_size)
  }
  else {
    plot(network, layout = coords, vertex.label.color = "black", edge.width = 3,
         vertex.label.dist = node_label_dist, vertex.size = vertex_size,
         vertex.label = node_labels)
  }
  title(main = title)
}

n <- 120
p <- 100

path <- paste0("Data/n", n, "_p", p, "/")

sim_random <- readRDS(file = paste0(path, "bdgraph_random_n", n, "_p", p, ".Rds"))
sim_bdgraph_sf <- readRDS(file = paste0(path, "bdgraph_scale-free_n", n, "_p", p, ".Rds"))
sim_huge_sf <- readRDS(file = paste0(path, "huge_scale-free_n", n, "_p", p, ".Rds"))
sim_hubs <- readRDS(file = paste0(path, "huge_hubs_n", n, "_p", p, ".Rds"))


setEPS()
postscript("Figures/Main_article/Simulated_structures.eps", width = 20, height = 5)

par(mfrow = c(1,4))
plot_network(sim_random[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_bdgraph_sf[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_huge_sf[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)
plot_network(sim_hubs[[1]]$theta, margins = c(0, 0, 0, 0), node_labels = NA)

dev.off()