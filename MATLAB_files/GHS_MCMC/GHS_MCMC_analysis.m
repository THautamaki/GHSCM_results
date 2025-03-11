clearvars;

% Initialize number of datasets.
n_datasets = 50;

% Initialize tables for all results and total times.
all_scores = table('Size', [n_datasets, 21], 'VariableTypes', ["double", "double", "double", "double", "double", "double", "double", ...
      "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"], ...
      'VariableNames', ["ACC", "ACC_bal", "MCC", "F1", "TPR", "TNR", "PPV", ...
      "NPV", "FPR", "FNR", "FDR", "FOR", "LRp", "LRn", "sl_omega", "f_norm", "f_norm_rel", "edge_count", "n", "round", "time"]);

total_times = table('Size', [0, 3], 'VariableTypes', ["double", "string", "double"], 'VariableNames', ["p", "structure", "time"]);

%% Initialize parfor.
y = ones(1,100);
parfor (i = 1:100)
     y(i) = i;
end
clear y;

%% Initialize needed parameters.
n = 120;
ps = [100, 200];
structures = ["random", "bdgraph_sf", "huge_sf", "hubs"];
path_tot_times = '..\..\Results_files\';

%% Run all analyzes.
for p = ps
    path_data = ['..\..\Data\n', num2str(n) '_p', num2str(p), '\'];
    path_res = ['..\..\Results_files\p' num2str(p), '\GHS_MCMC\'];
    for structure = structures
        t_start = tic;
        struct_char = char(structure);
        parfor r = 1:n_datasets
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
            [GHS_omega, ~] = GHS(S, n, 1000, 5000, 0);
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
            scores = table('Size', [1,21], 'VariableTypes', ["double", "double", "double", "double", "double", "double", "double", ...
                "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"], ...
                'VariableNames', ["ACC", "ACC_bal", "MCC", "F1", "TPR", "TNR", "PPV", ...
                "NPV", "FPR", "FNR", "FDR", "FOR", "LRp", "LRn", "sl_omega", "f_norm", "f_norm_rel", "edge_count", "n", "round", "time"]);
        
            scores(1, 1:14) = calculate_scores(cm);
            scores(1, "n") = {n};
            scores(1, "round") = {r};
            scores(1, "edge_count") = {cm{1,1} + cm{2,1}};
            sl_omega = stein_loss(omega, Omega_mean);
            scores(1, "sl_omega") = {sl_omega};
            scores(1, "f_norm") = {norm(omega - Omega_mean, "fro")};
            scores(1, "f_norm_rel") = {norm(omega - Omega_mean, "fro") / norm(omega, "fro")};
            scores(1, "time") = {time}
            all_scores(r,:) = scores;
            fprintf('Dataset %d finished \n',datanro);
        end
        total_times(end + 1, :) = {p, structure, toc(t_start)};
        writetable(all_scores, [path_res, struct_char, '_p', num2str(p), '_scores.csv'])
    end
end

%% Write total times into csv-file.
writetable(total_times, [path_tot_times, 'GHS_MCMC_total_times.csv'])


