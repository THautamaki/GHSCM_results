function cm = conf_matrix(truth, estimation)
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