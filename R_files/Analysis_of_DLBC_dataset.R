library(GHSGEM)

source("R_files/network_visualisation_and_scores.R")

######
# Combine the pathways, create a dataset with unique variables, transform it using huge.npn and save
# it as a csv file if data file not yet created.
if(!file.exists("Data/DLBC_dataset/rppadat_DLBC_npn.csv")) {
  load("Data/DLBC_dataset/rppadat_DLBC.rda")
  
  # Combine pathways.
  rppa_data <- data.frame(row.names = c(1:33))
  for (i in 1:12) {
    rppa_data <- cbind(rppa_data, rppadat[[i]])
  }
  
  # Select unique variable names.
  var_names <- unique(colnames(rppa_data))
  
  # Check that length of the variable names is 67.
  length(var_names)
  
  # Select variables and convert data to a matrix.
  rppa_data <- as.matrix(rppa_data[, var_names])
  
  # Data dimensions
  n <- nrow(rppa_data)
  p <- ncol(rppa_data)
  
  # Transform data using huge.npn().
  data_npn <- huge::huge.npn(rppa_data)
  
  # Write dataset to the csv-file.
  write.csv(data_npn, "Data/DLBC_dataset/rppadat_DLBC_npn.csv", row.names = FALSE)
}

######
# Read new, transformed dataset
data_npn <- as.matrix(read.csv(file = "Data/DLBC_dataset/rppadat_DLBC_npn.csv"))

# Data dimensions
n <- nrow(data_npn)
p <- ncol(data_npn)

# Run GHS GEM algorithm
GHSGEM_MAP <- GHS_MAP_estimation(data_npn, verbose = 1)

######
# Read adjacency matrices of other methods
GHS_MCMC_Theta <- as.matrix(read.csv("Results_files/DLBC_dataset/DLBC_data_npn_GHS_MCMC_Theta.txt",
                                     header = FALSE))

GHS_LLA_Theta <- as.matrix(read.csv("Results_files/DLBC_dataset/DLBC_data_npn_GHS_LLA_Theta.txt",
                                    header = FALSE))

GHSl_ECM_Theta <- as.matrix(read.csv("Results_files/DLBC_dataset/DLBC_data_npn_GHSl_ECM_Theta.txt",
                                     header = FALSE))

colnames(GHS_MCMC_Theta) <- rownames(GHS_MCMC_Theta) <- 1:p
colnames(GHS_LLA_Theta) <- rownames(GHS_LLA_Theta) <- 1:p
colnames(GHSl_ECM_Theta) <- rownames(GHSl_ECM_Theta) <- 1:p

######
# Calculate number of connections
sum(GHSGEM_MAP$Theta_est) / 2
sum(GHS_MCMC_Theta) / 2
sum(GHS_LLA_Theta) / 2
sum(GHSl_ECM_Theta) / 2

# Calculate node degrees
GHS_GEM_degrees <- colSums(GHSGEM_MAP$Theta_est)
GHS_MCMC_degrees <- colSums(GHS_MCMC_Theta)
GHS_LLA_degrees <- colSums(GHS_LLA_Theta)
GHSl_ECM_degrees <- colSums(GHSl_ECM_Theta)

# Calculate number of non-isolated nodes
sum(GHS_GEM_degrees > 0)
sum(GHS_MCMC_degrees > 0)
sum(GHS_LLA_degrees > 0)
sum(GHSl_ECM_degrees > 0)

######
# Calculate coordinates for the nodes
set.seed(13)
ceu_coords <- igraph::layout_with_fr(igraph::graph_from_adjacency_matrix(GHSGEM_MAP$Theta_est,
                                                                         mode = "undirected",
                                                                         diag = FALSE))

# Plot network estimates
setEPS()
postscript("Figures/Main_article/DLBC_network_estimates.eps", width = 20, height = 5)
par(mfrow = c(1,4))
plot_network(GHSGEM_MAP$Theta, layout = ceu_coords, node_size = GHS_GEM_degrees*1.5,
             node_labels = NA, margins = c(0,0,0,0))

plot_network(GHS_MCMC_Theta, layout = ceu_coords, node_size = GHS_MCMC_degrees*1.5,
             node_labels = NA, margins = c(0,0,0,0))

plot_network(GHS_LLA_Theta, layout = ceu_coords, node_size = GHS_LLA_degrees/1.5, node_labels = NA,
             margins = c(0,0,0,0))

plot_network(GHSl_ECM_Theta, layout = ceu_coords, node_size = GHSl_ECM_degrees/1.5, node_labels = NA,
             margins = c(0,0,0,0))
dev.off()

# Confusion matrix between GHS GEM and GHS MCMC estimates
conf_matrix(GHSGEM_MAP$Theta, GHS_MCMC_Theta)

# Set names for the nodes.
names(GHS_GEM_degrees) <- names(GHS_MCMC_degrees) <- 1:p

# Calculate how many nodes with the highest degree are the same
same_nodes <- c()
for (n1 in names(GHS_GEM_degrees[order(GHS_GEM_degrees, decreasing = TRUE)][1:10])) {
  for (n2 in names(GHS_MCMC_degrees[order(GHS_MCMC_degrees, decreasing = TRUE)][1:10])) {
    if (n1 == n2) {
      same_nodes <- c(same_nodes, n1)
    }
  }
}
same_nodes

# Highest degrees
GHS_GEM_degrees[order(GHS_GEM_degrees, decreasing = TRUE)][1:10]
GHS_MCMC_degrees[order(GHS_MCMC_degrees, decreasing = TRUE)][1:10]
GHS_LLA_degrees[order(GHS_LLA_degrees, decreasing = TRUE)][1:10]
GHSl_ECM_degrees[order(GHSl_ECM_degrees, decreasing = TRUE)][1:10]

