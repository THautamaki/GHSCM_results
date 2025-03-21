% Author: Tuomas Hautamäki, University of Oulu
function scores = calculate_scores(cm)
  % This function calculates performance scores for network estimation. The scores are same as binary
  % classification scores as the adjacency matrices are binary matrices.
  % 
  % Inputs
  %   cm  The 2 by 2 confusion matrix calculated using the function conf_matrix() or
  %       manually created 2 by 2 matrix using following order:
  %         TP  FN
  %         FP  TN
  %
  % Outputs:
  %   scores    A table, which contains following scores:
  %     ACC       Accuracy, which is (TP + TN) / (TP + TN + FN + FP).
  %     ACC_bal   Balanced accuracy, which is (TPR + TNR) / 2.
  %     MCC       Matthews correlation coefficient, which is (TP * TN - FP * FN) / sqrt{(TP + FP)(TP + FP)(TN + FP)(TN + FN)}.
  %     F1        F1-score, which is 2(PPV * TPR) / (PPV + TPR).
  %     TPR       True positive rate (or recall or sensitivity), which is TP / (TP + TN).
  %     TNR       True negative rate (or specificity or selectivity), which is TN / (TN + FP).
  %     PPV       Positive predictive value (or precision), which is TP / (TP + FP).
  %     NPV       Negative predictive value, which is TN / (TN + FN).
  %     FPR       False positive rate (type I error), which is 1 - TNR.
  %     FNR       False negative rate (type II error), which is 1 - TPR.
  %     FDR       False discovery rate, which is 1 - PPV.
  %     FOR       False omission rate, which is 1 - NPV.
  %     LRp       Positive likelihood ratio, which is TPR / FPR.
  %     LRn       Negative likelihood ratio, which is FNR / TNR.
  tp = cm{1,1};
  tn = cm{2,2};
  fp = cm{2,1};
  fn = cm{1,2};
  tpr = tp / (tp + fn);
  tnr = tn / (tn + fp);
  ppv = tp / (tp + fp);
  npv = tn / (tn + fn);
  fnr = 1 - tpr;
  fpr = 1 - tnr;
  fdr = 1 - ppv;
  FOR = 1 - npv;
  lr_plus = tpr / fpr;
  lr_neg = fnr / tnr;
  pt = sqrt(fpr) / (sqrt(tpr) + sqrt(fpr));
  ts = tp / (tp + fn + fp);
  acc = (tp + tn) / (tp + tn + fn + fp);
  bal_acc = (tpr + tnr) / 2;
  F1_score = 2 * (ppv * tpr) / (ppv + tpr);
  mcc = (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
  scores = table('Size', [1,14], 'VariableTypes', ["double", "double", "double", "double", "double", "double", "double", ...
    "double", "double", "double", "double", "double", "double", "double",], ...
    'VariableNames', ["ACC", "ACC_bal", "MCC", "F1", "TPR", "TNR", "PPV", ...
    "NPV", "FPR", "FNR", "FDR", "FOR", "LRp", "LRn"]);
  scores(1,1) = {acc};
  scores(1,2) = {bal_acc};
  scores(1,3) = {mcc};
  scores(1,4) = {F1_score};
  scores(1,5) = {tpr};
  scores(1,6) = {tnr};
  scores(1,7) = {ppv};
  scores(1,8) = {npv};
  scores(1,9) = {fpr};
  scores(1,10) = {fnr};
  scores(1,11) = {fdr};
  scores(1,12) = {FOR};
  scores(1,13) = {lr_plus};
  scores(1,14) = {lr_neg};
end