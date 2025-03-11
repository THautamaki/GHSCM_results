function loss = stein_loss(true_matrix, est_matrix)
  p = size(true_matrix, 1);
  P = est_matrix / true_matrix;
  loss = trace(P) - log(det(P)) - p;
end