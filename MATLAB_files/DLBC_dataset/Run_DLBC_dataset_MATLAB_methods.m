%% Read dataset.
clearvars;
rng("default") % Reset random number generator to default session setting.

datafile = ['..\..\Data\DLBC_dataset\rppadat_DLBC_npn.csv'];

data = readtable(datafile, "NumHeaderLines", 1);
data = table2array(data);

%% Extract dimensions and calculate a scatter matrix.
[n, p] = size(data);
S = data'*data;

%% Run MCMC algorithm.
[GHS_omega, ~] = GHS(S, n, 1000, 5000, 1);

%% Construct adjacency matrix using 50 percent credible intervals.
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

%% Write the adjacency matrix into text file.
writematrix(a_mat, "..\..\Results_files\DLBC_dataset\DLBC_data_npn_GHS_MCMC_Theta.txt")

%% Run GHS-like ECM. Initialize parameters.
[n1,q1] = size(data);
n_EMs = 50;
b_init = 0.02856817;
Omega_saves = zeros(p,p,n_EMs);

%% Generate initial values.
rng(20250326);
for i = 1:n_EMs % Do not use parfor here if results need to be reproducible.
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
    fprintf("Finished %d th start point generation out of %d start points \n", i, n_EMs);
    Omega_saves(:,:,i) = start_point;
end

%% Run ECM algorithm.
[Omega_est, ~] = Multi_start_point_Fixed_b_EM_HS_like(Omega_saves, S, n1, q1, b_init, n_EMs);

%% Calculate mean of the precision matrix.
Omega_mean = mean(Omega_est, 3);

%% Construct adjacency matrix.
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

%% Write the adjacency matrix into text file.
writematrix(a_mat, "..\..\Results_files\DLBC_dataset\DLBC_data_npn_GHSl_ECM_Theta.txt")

%% Run GHS LLA (Laplace). Find optimal value for tau.
desired_tau_ll = CV_HS_LLA_Laplace(data, 1);

%% Run LLA algorithm.
[Omega_est, ~] = Multi_start_point_Fixed_tau_HS_LLA_Laplace(data, desired_tau_ll, n_EMs, 20250327, 1);

%% Calculate mean of the precision matrix.
Omega_mean = mean(Omega_est, 3);

%% Construct adjacency matrix.
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

%% Write the adjacency matrix into text file.
writematrix(a_mat, "..\..\Results_files\DLBC_dataset\DLBC_data_npn_GHS_LLA_Theta.txt")

