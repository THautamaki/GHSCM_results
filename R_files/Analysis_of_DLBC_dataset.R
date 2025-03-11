library(huge)
library(igraph)
library(GHSGEM)

source("./network_visualisation_and_scores.R")

######
# Uncomment lines 12 to 37 if the DLBC dataset is not yet created.

# Combine the pathways, create a dataset with unique variables and save it as a csv file.
# Load rda-file.
# load("Data/DLBC_dataset/rppadat_DLBC.rda")
# 
# # Combine pathways.
# rppa_data <- data.frame(row.names = c(1:33))
# for (i in 1:12) {
#   rppa_data <- cbind(rppa_data, rppadat[[i]])
# }
# 
# # Select unique variable names.
# var_names <- unique(colnames(rppa_data))
# 
# # Check that length of the variable names is 67.
# length(var_names)
# 
# # Select variables and convert data to a matrix.
# rppa_data <- as.matrix(rppa_data[, var_names])
# 
# # Data dimensions
# n <- nrow(rppa_data)
# p <- ncol(rppa_data)
# 
# # Transform data using huge.npn().
# data_npn <- huge::huge.npn(rppa_data)
# 
# # Write dataset to the csv-file.
# write.csv(data_npn, "Data/DLBC_dataset/rppadat_DLBC_npn.csv", row.names = FALSE)

######
# Read new, transformed dataset
data_npn <- as.matrix(read.csv(file = "Data/CEU_parents_npn.csv"))

# Data dimensions
n <- nrow(rppa_data)
p <- ncol(rppa_data)

# Run GHS GEM algorithm
GHSGEM_MAP <- GHS_MAP_estimation(data_npn, verbose = 1)

######
# Read adjacency matrices of other methods
GHS_MCMC_Theta <- as.matrix(read.csv("Results_files/DLBC_dataset/DLBC_data_npn_GHS_MCMC_Theta_50_CI.txt",
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
sum(GHSGEM_MAP$Theta) / 2
sum(GHS_MCMC_Theta) / 2
sum(GHS_LLA_Theta) / 2
sum(GHSl_ECM_Theta) / 2

# Calculate node degrees
GHS_GEM_degrees <- colSums(GHSGEM_MAP$Theta)
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
ceu_coords <- igraph::layout_with_fr(igraph::graph_from_adjacency_matrix(GHSGEM_MAP$Theta,
                                                                         mode = "undirected",
                                                                         diag = FALSE))

# Plot network estimates
setEPS()
postscript("CEU_network_estimates.eps", width = 20, height = 5)
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




########

ghs_mcmc_theta <- as.matrix(read.csv("C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\Tuloksia\\DLBC_data_GHS_MCMC_Theta_50_CI.txt",
                                     header = FALSE))
ghs_mcmc_theta <- as.matrix(read.csv("C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\Tuloksia\\DLBC_data_scaled_GHS_MCMC_Theta_50_CI.txt",
                                     header = FALSE))
ghsl_ecm_theta <- as.matrix(read.csv("C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\Tuloksia\\DLBC_npn_data_GHSl_ECM_Theta.txt",
                                     header = FALSE))
ghs_lla_ll_theta <- as.matrix(read.csv("C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\Tuloksia\\DLBC_npn_data_GHS_LLA_ll_Theta.txt",
                                       header = FALSE))
ghs_lla_c_theta <- as.matrix(read.csv("C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\Tuloksia\\DLBC_npn_data_GHS_LLA_c_Theta.txt",
                                       header = FALSE))

gem_degrees <- colSums(bcell_map3$Theta)
#gem_degrees <- colSums(bcell_map2$Theta)
mcmc_degrees <- colSums(ghs_mcmc_theta)
ecm_degrees <- colSums(ghsl_ecm_theta)
lla_degrees <- colSums(ghs_lla_ll_theta)

names(gem_degrees) <- names(mcmc_degrees) <- 1:p

sum(bcell_map3$Theta) / 2
sum(bcell_map2$Theta) / 2
sum(ghs_mcmc_theta) / 2
sum(ghs_lla_c_theta) / 2
sum(ghs_lla_ll_theta) / 2
sum(ghsl_ecm_theta) / 2


set.seed(13)
net_coords <- igraph::layout_with_fr(graph_from_adjacency_matrix(bcell_map3$Theta, mode = "undirected", diag = F))

setEPS()
postscript("DLBC_network_estimates.eps", width = 20, height = 5)
par(mfrow = c(1,4))
plot_network(bcell_map3$Theta, layout = net_coords, margins = c(0, 0, 0, 0), node_labels = NA, node_size = gem_degrees, delete_isolates = FALSE)
plot_network(ghs_mcmc_theta, layout = net_coords, margins = c(0, 0, 0, 0), node_labels = NA, node_size = mcmc_degrees, delete_isolates = FALSE)
plot_network(ghs_lla_ll_theta, layout = net_coords, margins = c(0, 0, 0, 0), node_labels = NA, node_size = lla_degrees, delete_isolates = FALSE)
plot_network(ghsl_ecm_theta, layout = net_coords, margins = c(0, 0, 0, 0), node_labels = NA, node_size = ecm_degrees, delete_isolates = FALSE)
dev.off()


par(mfrow = c(1,2))
par(mar = c(5.1, 4.1, 2.1, 2.1))
hist(gem_degrees)
hist(mcmc_degrees)
hist(ecm_degrees)
hist(lla_degrees)

gem_degrees[order(gem_degrees, decreasing = TRUE)][1:10]
mcmc_degrees[order(mcmc_degrees, decreasing = TRUE)][1:10]
ecm_degrees[order(ecm_degrees, decreasing = TRUE)][1:10]
lla_degrees[order(lla_degrees, decreasing = TRUE)][1:10]

same <- c()
for (dg1 in names(gem_degrees[order(gem_degrees, decreasing = TRUE)][1:10])) {
  for (dg2 in names(mcmc_degrees[order(mcmc_degrees, decreasing = TRUE)][1:10])) {
    if (dg1 == dg2) {
      same <- c(same, dg1)
    }
  }
}

same

##########

conf_matrix(bcell_map3$Theta, ghs_mcmc_theta)

conf_matrix(bcell_map3$Theta, ghs_lla_c_theta)
conf_matrix(bcell_map3$Theta, ghs_lla_ll_theta)

conf_matrix(bcell_map3$Theta, ghsl_ecm_theta)

conf_matrix(ghs_lla_c_theta, ghsl_ecm_theta)
conf_matrix(ghs_lla_ll_theta, ghsl_ecm_theta)

#######
similarity_plots(bcell_map3$Theta, ghs_mcmc_theta)
