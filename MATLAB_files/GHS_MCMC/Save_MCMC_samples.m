clearvars;

% Initialize number of datasets.
n_datasets = 16;

runtimes = table('Size', [0, 4], 'VariableTypes', ["string", "double", "double", "double"], 'VariableNames', ["structure", "p", "datanro", "time"]);

%% Initialize parfor.
y = ones(1,100);
parfor i = 1:100
     y(i) = i;
end
clear y;

%% Initialize needed parameters.
n = 120;
ps = [100, 200];
structures = ["random", "bdgraph_sf", "huge_sf", "hubs"];
path_tot_times = '..\..\Results_files\';

datasets = [10 13 42 48 40  7 33  4 45 25  6 26  9 47 37 44];

%% Run all analyzes.
for p = ps
    path_data = ['..\..\Data\n', num2str(n) '_p', num2str(p), '\'];
    path_res = ['..\..\Results_files\p' num2str(p), '\GHS_MCMC\'];
    for structure = structures
        t_start = tic;
        struct_char = char(structure);
        parfor (r = 1:n_datasets, 16)
            datanro = datasets(r);
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
            [GHS_omega, ~] = GHS(S, n, 0, 10000, 0);
            time = toc(t2)
            
            runtimes(r, :) = {structure, p, datanro, time};
            s = struct("Omega", GHS_omega);
            save([path_res, struct_char, '\', struct_char, '_p', num2str(p), '_Omega_samples_datanro_', num2str(datanro), '.mat'], "-fromstruct", s)
            fprintf('Dataset %d finished \n',datanro);
        end
        writetable(runtimes, [path_res, struct_char, '_p', num2str(p), '_runtimes.csv'])
    end
end
