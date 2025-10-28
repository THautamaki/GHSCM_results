p <- 2000
n <- 0.6 * p

sim <- huge::huge.generator(n = n, d = p, graph = "scale-free")

map <- GHSCM::GHS_MAP_estimation(sim$data, max_iterations = 2, verbose = 2)