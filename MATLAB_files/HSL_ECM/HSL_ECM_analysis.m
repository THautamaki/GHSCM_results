clearvars;

% Initialize number of datasets and number of initial values.
n_datasets = 50;
n_EMs = 50;

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

%%
for p = ps
    path_data = ['..\..\Data\n', num2str(n) '_p', num2str(p), '\'];
    path_res = ['..\..\Results_files\p' num2str(p), '\HSL_ECM\'];
    Omega_saves = zeros(p,p,n_EMs);
    parfor i = 1:n_EMs % choose parfor instead of "for" loop if running on a multi core CPU
        start_point = eye(p);
        for row = 2:p
            d = 0;
                while d~=p
                    row_seq = row:1:p;
                    col_seq = 1:1:(p-row+1);
    
                    rand_noise = -0.05 + rand(1, length(row_seq))*2*0.05;
                    
                    lin_idcs = sub2ind(size(start_point), row_seq, col_seq);
                    start_point(lin_idcs) = rand_noise;
    
                    lin_idcs = sub2ind(size(start_point), col_seq, row_seq);
                    start_point(lin_idcs) = rand_noise;
    
                    d = eig(start_point);
                    d = sum(d>0);
                end
        end
        Omega_saves(:,:,i) = start_point;
    end
    if(p == 100 && n == 120)
        fixed_b = 0.0143;
    elseif(p == 200 && n == 120)
        fixed_b = 0.0169;
    end
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
            
            S = data'*data;

            t3 = tic;
            [Omega_est, ~, ~, ~] = Multi_start_point_Fixed_b_EM_HS_like(Omega_saves, S, n, p, fixed_b, 50);
            time = toc(t3);

            Omega_mean = mean(Omega_est, 3);
            
            a_mat = zeros(p);
            for i = 1:p
                for j = i:p
                    if i == j
                        continue
                    elseif Omega_mean(i,j) == 0
                        continue
                    else
                        a_mat(i,j) = 1;
                        a_mat(j,i) = 1;
                    end
                end
            end
            
            cm = conf_matrix(theta, a_mat);
            scores = table('Size', [1,21], 'VariableTypes', ["double", "double", "double", "double", "double", "double", "double", "double", ...
                "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"], ...
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
            scores(1, "time") = {time};
            all_scores(r,:) = scores;
            fprintf('Dataset %d finished \n',datanro);
        end
        total_times(end + 1, :) = {p, structure, toc(t_start)};
        writetable(all_scores, [path_res, struct_char, '_p', num2str(p), '_scores.csv'])
    end
end

%% Write total times into csv-file.
writetable(total_times, [path_tot_times, 'HSL_ECM_total_times.csv'])


