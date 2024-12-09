clearvars;

total_sets = 4;

all_scores = table('Size', [total_sets,20], 'VariableTypes', ["double", "double", "double", "double", "double", "double", "double", ...
      "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"], ...
      'VariableNames', ["ACC", "ACC_bal", "MCC", "F1", "TPR", "TNR", "PPV", ...
      "NPV", "FPR", "FNR", "FDR", "FOR", "LRp", "LRn", "sl_omega", "f_norm", "edge_count", "n", "round", "time"]);

%%
n = 120;
%ps = [100, 200];
ps = 100;
%structures = ["random", "bdgraph_sf", "huge_sf", "hubs"];
structures = ["random", "bdgraph_sf"];

%%
total_times = table('Size', [0, 3], 'VariableTypes', ["double", "string", "double"], 'VariableNames', ["p", "structure", "time"]);
for p = ps
    %path_data = ['C:\Users\tume8\Documents\data\Artikkeli\p', num2str(p), '\'];
    path_data = ['C:\Users\thautama\OneDrive - Oulun yliopisto\Documents\data\Artikkeli\p', num2str(p), '\'];
    %path_res = ['C:\Users\tume8\Documents\Tuloksia\Artikkeli\p' num2str(p), '\GHS_MCMC\'];
    path_res = ['C:\Users\thautama\OneDrive - Oulun yliopisto\Documents\Tuloksia\Artikkeli\p' num2str(p), '\GHS_MCMC\'];
    for structure = structures
        t_start = tic;
        struct_char = char(structure);
        parfor r = 1:total_sets
            datanro = r;
            fprintf('Dataset %d is in process.\n',datanro);
            datafile = [path_data, struct_char, '\', struct_char, '_data_nro_', num2str(datanro), '.csv'];
            covfile = [path_data, struct_char, '\', struct_char, '_sigma_nro_', num2str(datanro), '.csv'];
            precfile = [path_data, struct_char, '\', struct_char, '_omega_nro_', num2str(datanro), '.csv'];
            adjfile = [path_data, struct_char, '\', struct_char, '_theta_nro_', num2str(datanro), '.csv'];
        
            data = readtable(datafile, "VariableNamesRow", 1);
            data = table2array(data);
            [n1,q1] = size(data);
        
            sigma = readtable(covfile,  "VariableNamesRow", 1);
            sigma = table2array(sigma);
            
            omega = readtable(precfile, "VariableNamesRow", 1);
            omega = table2array(omega);
            
            theta = readtable(adjfile,  "VariableNamesRow", 1);
            theta = table2array(theta);
        
            S = data'*data;
        
            t2 = tic;
            [GHS_omega, ~] = GHS(S, n, 100, 500, 0);
            time = toc(t2)
            Omega_mean = mean(GHS_omega, 3);
            
            quant = 0.25;
            low_q = quantile(GHS_omega, quant, 3);
            up_q = quantile(GHS_omega, 1 - quant, 3);
            
            a_mat = zeros(p);
            for i = 1:p
                for j = i:p
                    if i == j
                        continue
                    elseif low_q(i,j) <= 0 && up_q(i,j) > 0
                        continue
                    else
                        a_mat(i,j) = 1;
                        a_mat(j,i) = 1;
                    end
                end
            end
        
            cm = conf_matrix(theta, a_mat);
            scores = table('Size', [1,20], 'VariableTypes', ["double", "double", "double", "double", "double", "double", "double", ...
                "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"], ...
                'VariableNames', ["ACC", "ACC_bal", "MCC", "F1", "TPR", "TNR", "PPV", ...
                "NPV", "FPR", "FNR", "FDR", "FOR", "LRp", "LRn", "sl_omega", "f_norm", "edge_count", "n", "round", "time"]);
        
            scores(1, 1:14) = calculate_scores(cm);
            scores(1, "n") = {n};
            scores(1, "round") = {r};
            scores(1, "edge_count") = {cm{1,1} + cm{2,1}};
            sl_omega = stein_loss(omega, Omega_mean);
            scores(1, "sl_omega") = {sl_omega};
            scores(1, "f_norm") = {norm(omega - Omega_mean, "fro")};
            scores(1, "time") = {time}
            all_scores(r,:) = scores;
            fprintf('Dataset %d finished \n',datanro);
        end
        total_times(end + 1, :) = {p, structure, toc(t_start)};
        writetable(all_scores, [path_res, struct_char, '_p', num2str(p), '_scores.csv'])
    end
end



