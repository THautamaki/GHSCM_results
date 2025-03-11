library(igraph)
library(viridis)
library(plot.matrix)

stein_loss <- function(m1, m2){
  n <- nrow(m1)
  P <- chol2inv(chol(m2)) %*% m1
  return(sum(diag(P)) - log(det(P)) - n)
}

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

conf_matrix <- function(truth, estimation, margins = FALSE, normalize = FALSE,
                        undirected = TRUE) {
  same_edges <- truth * estimation
  diff <- truth - estimation
  summ <- truth + estimation
  p <- dim(truth)[1]
  max_edges <- (p^2 - p)
  tp <- sum(same_edges)
  tn <- sum(summ == 0) - p
  fp <- sum(diff == -1)
  fn <- sum((same_edges - truth) == -1)
  P <- sum(truth)
  N <- max_edges - P
  EP <- sum(estimation)
  EN <- max_edges - EP
  cm <- matrix(c(tp, fn, fp, tn), nrow = 2, byrow = TRUE,
               dimnames = list(c("True P", "True N"), c("Estim. P", "Estim. N")))
  if (margins) {
    cm <- matrix(c(tp, fn, P, fp, tn, N, EP, EN, max_edges), nrow = 3, byrow = TRUE,
                 dimnames = list(c("True P", "True N", "Sum"), c("Estim. P", "Estim. N", "Sum")))
  }
  if (undirected) {
    cm <- cm * 0.5
  }
  if (normalize) {
    cm <- matrix(c(tp/P, fn/P, fp/N, tn/N), nrow = 2, byrow = TRUE,
                 dimnames = list(c("True P", "True N"), c("Estim. P", "Estim. N")))
  }
  return(cm)
}

calculate_scores <- function(cm) {
  tp <- cm[1,1]
  tn <- cm[2,2]
  fp <- cm[2,1]
  fn <- cm[1,2]
  tpr <- tp / (tp + fn)
  tnr <- tn / (tn + fp)
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)
  fnr <- 1 - tpr
  fpr <- 1 - tnr
  fdr <- 1 - ppv
  FOR <- 1 - npv
  lr_plus <- tpr / fpr
  lr_neg <- fnr / tnr
  pt <- sqrt(fpr) / (sqrt(tpr) + sqrt(fpr))
  ts <- tp / (tp + fn + fp)
  fm <- sqrt(ppv * tpr)
  mk <- ppv + npv - 1
  acc <- (tp + tn) / (tp + tn + fn + fp)
  bal_acc <- (tpr + tnr) / 2
  F1_score <- 2 * (ppv * tpr) / (ppv + tpr)
  mcc <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  results <- data.frame(ACC = acc, ACC_bal = bal_acc, MCC = mcc, F1 = F1_score,
                        TPR = tpr, TNR = tnr, PPV = ppv, NPV = npv, FPR = fpr,
                        FNR = fnr, FDR = fdr, FOR = FOR, PT = pt, TS = ts, FM = fm, MK = mk,
                        LRp = lr_plus, LRn = lr_neg)
  return(results)
}


similarity_plots <- function(adjacency_T, adjacency_E, delete_isolates = FALSE) {
  par(mfrow = c(2,2))
  samat <- adjacency_T * adjacency_E
  plot_network(samat, "True positives", delete_isolates = delete_isolates)
  fn <- adjacency_T - samat 
  plot_network(fn, "False negatives", delete_isolates = delete_isolates)
  fp <- adjacency_E - samat 
  plot_network(fp, "False positives", delete_isolates = delete_isolates)
  plot_network(adjacency_T, "True network", delete_isolates = delete_isolates)
}

adjacency_matrix <- function(coef_matrix, condition = "AND") {
  p <- dim(coef_matrix)[1]
  adj_matrix <- matrix(0, nrow = p, ncol = p)
  if (condition != "AND" && condition != "OR") {
    print("Wrong condition! Only AND and OR are accepted.")
    return(adj_matrix)
  }
  for (i in 1:p) {
    for (j in i:p) {
      if (i == j) {
        next
      }
      if (condition == "AND") {
        if (coef_matrix[i, j] != 0 && coef_matrix[j, i] != 0) {
          adj_matrix[i, j] <- 1
          adj_matrix[j, i] <- 1
        }
      } 
      else if (condition == "OR") {
        if (coef_matrix[i, j] != 0 || coef_matrix[j, i] != 0) {
          adj_matrix[i, j] <- 1
          adj_matrix[j, i] <- 1
        }
      }
    }
  }
  return(adj_matrix)
}

plot_matrix <- function(mat, title = "", breaks = c()) {
  par(mar = c(4,4,2,3))
  key_params <- list(side = 4, font = 1, cex.axis = 0.8)
  if (length(breaks) == 0) {
    plot(mat, main = title, border = NA, col = magma, key = key_params, fmt.key = "%.2f",
         spacing.key = c(3, 2, 2))
  }
  else {
    plot(mat, main = title, border = NA, col = magma, key = key_params, fmt.key = "%.2f",
         spacing.key = c(3, 2, 2), breaks = breaks)
  }
}

partial_corr <- function(Omega) {
  p <- dim(Omega)[1]
  part_corr <- matrix(NA, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in i:p) {
      if (i != j)
        part_corr[i,j] <- part_corr[j,i] <- - Omega[i,j] / sqrt(Omega[i,i] * Omega[j,j]) 
    }
  }
  return(part_corr)
}
