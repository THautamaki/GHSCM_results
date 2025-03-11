library(huge)
library(igraph)
library(GHSGEM)

source("./network_visualisation_and_scores.R")

##############

load("C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\data\\TCGA_RPPA\\rppadat\\rppadat_DLBC.rda")

rppadat

rppadat[[1]]

rppa_data <- data.frame(row.names = c(1:33))
for (i in 1:12) {
  rppa_data <- cbind(rppa_data, rppadat[[i]])
}
rppa_data

var_names <- unique(colnames(rppa_data))

length(var_names)

#var_names == colnames(rppa_data)

rppa_data <- as.matrix(rppa_data[, var_names])

#rppa_data <- RPPAdat

n <- nrow(rppa_data)
p <- ncol(rppa_data)

data_norm <- huge::huge.npn(rppa_data)
#data_cent <- scale(rppa_data, center = TRUE, scale = FALSE)
#data_scale <- scale(rppa_data, center = TRUE, scale = TRUE)

par(mfrow = c(5,5))
par(mar = c(5.1, 4.1, 2.1, 2.1))
for (i in 51:67) {
  plot(density(data_norm[,i]))
}
for (i in 51:67) {
  plot(density(data_scale[,i]))
}

#write.csv(data_norm, "C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\data\\TCGA_RPPA\\rppadat_DLBC_npn.csv", row.names = F)
#write.csv(data_norm, "C:\\Users\\thautama\\OneDrive - Oulun yliopisto\\Documents\\data\\TCGA_RPPA\\rppadat_DLBC_scaled.csv", row.names = F)


?huge.npn

#data_cent <- scale(RPPAdat, center = TRUE, scale = TRUE)

#####

hist(apply(data_cent, 2, var))

par(mfrow = c(5,5))
for (i in 1:25) {
  plot(density(data_cent[,i]))
  lines(density(data_norm[,i]), col = "red")
  lines(density(data_scale[,i]), col = "blue")
}

#####

start <- Sys.time()
beam_rppa_est <- beam(data_norm)
beam_rppa_sel <- beam.select(beam_rppa_est, thres = 0.1)
end <- Sys.time()
end - start

beam_theta <- as.matrix(as_adjacency_matrix(ugraph(beam_rppa_sel), names = FALSE))

par(mfrow = c(1,1))
plot_network(beam_theta)

sum(beam_theta) / 2

colSums(beam_theta)

#########

bcell_map <- GHS_MAP_estimate(data_cent, verbose = 1, p0 = 100)

plot_network(bcell_map$Theta)

sum(bcell_map$Theta) / 2

colSums(bcell_map$Theta)

########

bcell_map2 <- GHS_MAP_estimate(data_scale, verbose = 1)

par(mfrow = c(1,1))
plot_network(bcell_map2$Theta)

sum(bcell_map2$Theta) / 2

hist(colSums(bcell_map2$Theta))

########

bcell_map3 <- GHS_MAP_estimation(data_norm, verbose = 1)

par(mfrow = c(1,1))
plot_network(bcell_map3$Theta)

sum(bcell_map3$Theta) / 2

hist(colSums(bcell_map3$Theta))


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
