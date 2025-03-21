% Author: Tuomas Hautamäki, University of Oulu
function cm = conf_matrix(truth, estimation)
  % Construct confusion matrix using the true adjacency and estimated adjacency matrices of
  % the network. Can be used both directed and undirected networks.
  %
  % Inputs
  %   truth       The true p by p adjacency matrix of the network.
  %   estimation  The estimated \code{p} by \code{p} adjacency matrix of the network.
  % 
  % Outputs
  %   cm          A 2 by 2 table, where the first element is number of the true positives (TP),
  %               second element is number of the false negatives (FN), third element is number of
  %               the falste positives (FP) and last element is number of the true negatives (TN).
  same_edges = truth .* estimation;
  diff = truth - estimation;
  summ = truth + estimation;
  p = size(truth, 1);
  max_edges = (p^2 - p);
  tp = sum(same_edges, "all");
  tn = sum(summ == 0, "all") - p;
  fp = sum(diff == -1, "all");
  fn = sum((same_edges - truth) == -1, "all");
  P = sum(truth, "all");
  N = max_edges - P;
  EP = sum(estimation, "all");
  EN = max_edges - EP;
  cm = table('Size', [2,2], 'VariableTypes', ["double", "double"], 'VariableNames', ["Estim. P", "Estim. N"], 'RowNames', ["True P", "True N"]);
  cm(1,1) = {tp};
  cm(1,2) = {fn};
  cm(2,1) = {fp};
  cm(2,2) = {tn};
  cm = cm .* 0.5;
end