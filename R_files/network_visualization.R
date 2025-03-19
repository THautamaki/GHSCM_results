library(igraph)

plot_network <- function(a_mat, title = "", layout = NULL, node_size = 3, node_label_dist = 0.8,
                         node_labels = NULL, delete_isolates = FALSE, delete_nodes_degree = 0, 
                         margins = c(5.1, 4.1, 4.1, 2.1)) {
  par(mar = margins)
  network <- graph_from_adjacency_matrix(a_mat, mode = "undirected", diag = F)
  if (delete_isolates) {
    network <- delete_vertices(network, V(network)[degree(network) == 0])
  }
  if (delete_nodes_degree > 0) {
    network <- delete_vertices(network, V(network)[degree(network) <= delete_nodes_degree])
  }
  if (is.null(layout)) {
    coords <- layout_with_fr(network)
  }
  else {
    coords <- layout
  }
  if (is.null(node_labels)) {
    plot(network, layout = coords, vertex.label.color = "black", edge.width = 3,
         vertex.label.dist = node_label_dist, vertex.size = node_size)
  }
  else {
    plot(network, layout = coords, vertex.label.color = "black", edge.width = 3,
         vertex.label.dist = node_label_dist, vertex.size = node_size,
         vertex.label = node_labels)
  }
  title(main = title)
}

similarity_plots <- function(adjacency_T, adjacency_E, delete_isolates = FALSE) {
  par(mfrow = c(2,2))
  same <- adjacency_T * adjacency_E
  plot_network(same, "True positives", delete_isolates = delete_isolates)
  fn <- adjacency_T - same 
  plot_network(fn, "False negatives", delete_isolates = delete_isolates)
  fp <- adjacency_E - same 
  plot_network(fp, "False positives", delete_isolates = delete_isolates)
  plot_network(adjacency_T, "True network", delete_isolates = delete_isolates)
}