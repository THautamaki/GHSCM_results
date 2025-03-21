% Author: Ksheera Sagar K. N., Purdue University
% Edits: Tuomas Hautamäki, University of Oulu
%	- Added input "b_init", which can be calculated using the R script scale_factor_computation.R
%	  from Sagar's GitHub page https://github.com/sagarknk/Graphical_HSL.

function [Omega_est, final_nu_matrix, total_iterations, each_time_taken]  = ...
    Multi_start_point_Fixed_b_EM_HS_like(Omega_saves, S, n, q, b_init, n_EMs)

    Omega_est = zeros(q,q,n_EMs);
    final_nu_matrix = zeros(q,q,n_EMs);
    total_iterations = zeros (1, n_EMs);
    each_time_taken = zeros(1, n_EMs);
    
   parfor init_iter = 1:n_EMs %choose parfor instead of "for" loop if running on a multi core CPU
       
    Omega_init = Omega_saves(:,:,init_iter); % Initialising the start of EM algorithm
   
    %%%%% Fixed scale paramter estimated using Piironen and Vehtari (2017)
    %if(p == 100 && n == 120)
    %    b_init = 0.0143;
    %elseif(p == 200 && n == 120)
    %    b_init = 0.0169;
    %end
    
    %%% ind_all is a matrix which contains q-1 rows and q columns.
    %%% k^th column of ind_all contains {1,2,...,q}\{k}.
    
    ind_all = zeros(q-1,q);
    for i = 1:q
           if i==1  
               ind = (2:q)'; 
          elseif i==q
              ind = (1:q-1)'; 
           else
               ind = [1:i-1,i+1:q]';
           end
           ind_all(:,i) = ind;
    end

    nu_matrix = zeros(q,q); % latent parameter matrix
    Omega_current  = Omega_init; % current value of precision matrix
    Omega_next = eye(q); % setting the updated value of precision matrix to Identity
    norm_diff = norm(Omega_current - Omega_next, 'fro'); % calculating
    %%% forbenious norm  between the current state and updated state
    
    iter = 1; % Starting the iterations with iter = 1
    time_start = tic; % start measuring the time 
    
    while norm_diff > 1e-3

        Omega_current = Omega_init;
        nu_matrix = zeros(q,q);
        
        
        for i = 2:q
            for j= 1:(i-1)
                nu_matrix(i,j) = (1/(log(1+(b_init/Omega_init(i,j)^2))))*((b_init)^(2)/(( Omega_init(i,j)^2 + b_init)*(Omega_init(i,j))^2));
            end 
        end 
        
        %%%% Numerical Stability %%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		temp_bool = isnan(nu_matrix);
        nu_matrix(temp_bool) = Inf;
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Next line is to make sure that both 
        %%% upper and lower triangular entries have same penalty
        nu_matrix = nu_matrix + nu_matrix';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         for i = 1:q

              ind = ind_all(:,i);

              Omega_11 = Omega_init(ind,ind); 
              s_12 = S(ind,i); 
              s_22 = S(i,i);
              nu_matrix_elements_col = nu_matrix(ind, i);
              D_matrix = 0.5*diag(b_init ./ nu_matrix_elements_col);

              gamma_hat = n/s_22;

              inv_Omega_11 = inv(Omega_11);

              try 
                  chol_s22_inv_Omega_11 = chol((s_22 * inv_Omega_11));
              catch 
                  break;
              end 

              inverse_of_A = D_matrix - D_matrix*chol_s22_inv_Omega_11' ...
                  *((chol_s22_inv_Omega_11*D_matrix*chol_s22_inv_Omega_11' + eye(q-1))\...
                  chol_s22_inv_Omega_11*D_matrix);  

              beta_hat = -inverse_of_A*s_12 ; 

              %%% update Omega
              Omega_init(i,ind) = beta_hat; 
              Omega_init(ind,i) = beta_hat;
              Omega_init(i,i) = gamma_hat + beta_hat'*inv_Omega_11* beta_hat;
         end
         
         Omega_next = Omega_init;
         norm_diff = norm(Omega_next - Omega_current, 'fro');
         iter = iter +1 ;
         
         % while loop of EM ends 
    end 
    
    each_time_taken(1,init_iter) = toc(time_start);
    final_nu_matrix(:,:,init_iter) = nu_matrix;
    total_iterations(1, init_iter) = iter-1;
    Omega_est(:,:,init_iter) = Omega_init;

   end 
end 