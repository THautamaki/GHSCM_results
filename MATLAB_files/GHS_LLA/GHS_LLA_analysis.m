clearvars;

% Initialize number of datasets.
n_datasets = 50;

% Initialize tables for all results and total times.
all_scores = table('Size', [n_datasets, 22], 'VariableTypes', ["double", "double", "double", "double", "double", "double", "double", "double", ...
      "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"], ...
      'VariableNames', ["ACC", "ACC_bal", "MCC", "F1", "TPR", "TNR", "PPV", ...
      "NPV", "FPR", "FNR", "FDR", "FOR", "LRp", "LRn", "sl_omega", "f_norm", "f_norm_rel", "edge_count", "n", "round", "time", "tau_time"]);

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
    path_res = ['..\..\Results_files\p' num2str(p), '\GHS_LLA\'];
    for structure = structures
        t_start = tic;
        struct_char = char(structure);
        for r = 1:n_datasets
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
            
            %S = data'*data;
            t2 = tic;
            desired_tau_ll = CV_HS_LLA_Laplace(data, 0);
            tau_calc_time = toc(t2);
            if length(desired_tau_ll) > 1
                desired_tau_ll = desired_tau_ll(1,end);
            end

            t3 = tic;
            [Omega_est_ll, ~, ~] = Multi_start_point_Fixed_tau_HS_LLA_Laplace(data, desired_tau_ll, 50, 0);
            time = toc(t3);

            Omega_mean_ll = mean(Omega_est_ll, 3);
            
            a_mat_ll = zeros(p);
            for i = 1:p
                for j = i:p
                    if i == j
                        continue
                    elseif Omega_mean_ll(i,j) == 0
                        continue
                    else
                        a_mat_ll(i,j) = 1;
                        a_mat_ll(j,i) = 1;
                    end
                end
            end
            
            cm = conf_matrix(theta, a_mat_ll);
            scores = table('Size', [1,22], 'VariableTypes', ["double", "double", "double", "double", "double", "double", "double", "double", ...
                "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"], ...
                'VariableNames', ["ACC", "ACC_bal", "MCC", "F1", "TPR", "TNR", "PPV", ...
                "NPV", "FPR", "FNR", "FDR", "FOR", "LRp", "LRn", "sl_omega", "f_norm", "f_norm_rel", "edge_count", "n", "round", "time", "tau_time"]);
        
            scores(1, 1:14) = calculate_scores(cm);
            scores(1, "n") = {n};
            scores(1, "round") = {r};
            scores(1, "edge_count") = {cm{1,1} + cm{2,1}};
            sl_omega = stein_loss(omega, Omega_mean_ll);
            scores(1, "sl_omega") = {sl_omega};
            scores(1, "f_norm") = {norm(omega - Omega_mean_ll, "fro")};
            scores(1, "f_norm_rel") = {norm(omega - Omega_mean_ll, "fro") / norm(omega, "fro")};
            scores(1, "time") = {time};
            scores(1, "tau_time") = {tau_calc_time};
            all_scores(r,:) = scores;
            fprintf('Dataset %d finished \n',datanro);
        end
        total_times(end + 1, :) = {p, structure, toc(t_start)};
        writetable(all_scores, [path_res, struct_char, '_p', num2str(p), '_scores.csv'])
    end
end

%% Write total times into csv-file.
writetable(total_times, [path_tot_times, 'GHS_LLA_ll_total_times.csv'])


