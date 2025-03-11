library(GHSGEM)

source("network_visualisation_and_scores.R")

######
# Create dataset with 100 variables and save it as csv file

# data1 <- read.csv("C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\data\\genevar\\CEU_parents-normalised.csv")
# 
# varnames <- data1[,1]
# data1 <- data1[,-1]
# data1 <- t(data1)
# colnames(data1) <- varnames
# 
# variances <- apply(data1, 2, var)
# sum(variances > 1.35)
# 
# data2 <- data1[, variances > 1.35]
# 
# data_npn <- huge.npn(data2)

######
# Read new, transformed dataset
data_npn <- as.matrix(read.csv(file = "C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\data\\genevar\\CEU_parents_npn.csv"))

# Data dimensions
n <- nrow(data_npn)
p <- ncol(data_npn)

######
# Run GHS GEM algorithm
GHSGEM_MAP <- GHS_MAP_estimation(data_npn, verbose = 1, max_iterations = 500)

######
# Read adjacency matrices of other methods
GHS_MCMC_Theta <- as.matrix(read.csv("C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\Tuloksia\\CEU_data_npn_GHS_MCMC_Theta.txt",
                                     header = FALSE))

GHS_LLA_Theta <- as.matrix(read.csv("C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\Tuloksia\\CEU_data_npn_GHS_LLA_Theta.txt",
                                    header = FALSE))

GHSl_ECM_Theta <- as.matrix(read.csv("C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\Tuloksia\\CEU_data_npn_GHSl_ECM_Theta.txt",
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
set.seed(8)
ceu_coords <- igraph::layout_with_fr(igraph::graph_from_adjacency_matrix(GHSGEM_MAP$Theta,
                                                                         mode = "undirected",
                                                                         diag = F))

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
for (n1 in names(GHS_GEM_degrees[order(GHS_GEM_degrees)][95:100])) {
  for (n2 in names(GHS_MCMC_degrees[order(GHS_MCMC_degrees)][95:100])) {
    if (n1 == n2) {
      same_nodes <- c(same_nodes, n1)
    }
  }
}
same_nodes

# Highest degrees
GHS_GEM_degrees[order(GHS_GEM_degrees)][95:100]
GHS_MCMC_degrees[order(GHS_MCMC_degrees)][95:100]
GHS_LLA_degrees[order(GHS_LLA_degrees)][95:100]
GHSl_ECM_degrees[order(GHSl_ECM_degrees)][95:100]
