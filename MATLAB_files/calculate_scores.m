function scores = calculate_scores(cm) 
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