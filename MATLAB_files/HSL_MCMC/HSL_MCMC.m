% Sampling for graphical horseshoe-like
% Original author: Ksheera Sagar K. N., Purdue University
% Edits: Tuomas Hautamäki, University of Oulu
%	- Removed T_sq_vector and tau_sq_save outputs.
%	- Added input "verbose".

function [omega_save,tau_sq_save] = HSL_MCMC(S,n,burnin,nmc,verbose)
% GHS-like MCMC sampler using data-augmented
% block (column-wise) Gibbs sampler
% Input:
%     S = Y'*Y:    sample covariance matrix * n
%     n:		   sample size
%     burnin, nmc: number of MCMC burnins and saved samples
%	  verbose;	   if 1, prints additional information every 100th iteration

% Output:
%     omega_save:  p by p by nmc matrices of saved posterior samples of
%                  precision matrix
%     omega_vector_save: Removed! p*(p-1)/2 by nmc matrices of saved posterior samples of
%                        lower(upper) trinagular entries of the precision matrix
%     T_sq_save:   Removed! p*(p-1)/2 by nmc vector of saved samples of t_{ij}
%                  squared (local tuning parameter)
%     tau_sq_save: 1 by nmc vector of saved samples of tau squared (global
%                  tuning parameter)

[p] = size(S,1); indmx = reshape([1:p^2],p,p); 
upperind = indmx(triu(indmx,1)>0); 

indmx_t = indmx';
lowerind = indmx_t(triu(indmx_t,1)>0); 

omega_save = zeros(p,p,nmc);
%omega_vector_save = zeros(p*(p-1)*0.5, nmc);
tau_sq_save = zeros(1, nmc);
%T_sq_save = zeros(p*(p-1)/2,nmc);

ind_noi_all = zeros(p-1,p);
for i = 1:p
       if i==1  
       ind_noi = [2:p]'; 
      elseif i==p
       ind_noi = [1:p-1]'; 
      else
       ind_noi = [1:i-1,i+1:p]';
       end
       
       ind_noi_all(:,i) = ind_noi;
end

% set initial values
Omega = eye(p); Sigma = eye(p);
T(1:p,1:p) = 1; M(1:p,1:p) = 1; tau_sq = 1; xi = 1;

%tic;
for iter = 1: burnin+nmc
	if(verbose > 0)
		if(mod(iter,100)==0)
			fprintf('iter = %d \n',iter);
		end
	end
	
%%% sample Sigma and Omega=inv(Sigma)
    for i = 1:p
      ind_noi = ind_noi_all(:,i);     
      Sigma_11 = Sigma(ind_noi,ind_noi); sigma_12 = Sigma(ind_noi,i);
      sigma_22 = Sigma(i,i);
      s_21 = S(ind_noi,i); s_22 = S(i,i);
      T_sq_12 = T(ind_noi,i); M_12 = M(ind_noi,i);
      %% sample gamma and beta
      gamma = gamrnd((n/2+1),2/s_22);    % random gamma with shape=n/2+1, rate=s_22/2
      inv_Omega_11 = Sigma_11 - sigma_12*sigma_12'/sigma_22;
      inv_C = s_22*inv_Omega_11+diag(T_sq_12./tau_sq);
      inv_C_chol = chol(inv_C);
      mu_i = -inv_C\s_21;
      beta = mu_i+ inv_C_chol\randn(p-1,1);
      omega_12 = beta; omega_22 = gamma + beta'*inv_Omega_11*beta;
      %% sample lambda_sq and nu
      rate = omega_12.^2/(2*tau_sq)+(M_12 ./ 2);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      T_sq_12_new = gamrnd(3/2,1./rate); 
      T_sq_12(T_sq_12_new<1e15)=T_sq_12_new(T_sq_12_new<1e15);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      M_12 = -1.*(2./T_sq_12).*log(1-rand([p-1,1]).*(1-exp(-0.5.*T_sq_12)));
      %% update Omega, Sigma, Lambda_sq, Nu
      Omega(i,ind_noi) = omega_12; Omega(ind_noi,i) = omega_12;
      Omega(i,i) = omega_22;
      temp = inv_Omega_11*beta;
      Sigma_11 = inv_Omega_11 + temp*temp'/gamma;
      sigma_12 = -temp/gamma; sigma_22 = 1/gamma;
      Sigma(ind_noi,ind_noi) = Sigma_11; Sigma(i,i) = sigma_22;
      Sigma(i,ind_noi) = sigma_12; Sigma(ind_noi,i) = sigma_12;
      T(i,ind_noi) = T_sq_12; T(ind_noi,i) = T_sq_12;
      M(i,ind_noi) = M_12; M(ind_noi,i) = M_12;
    end
    
%%% sample tau_sq and xi
    omega_vector = Omega(tril(true(size(Omega)),-1));
    T_sq_vector = T(tril(true(size(T)),-1));
    rate = 1/xi + sum((omega_vector.^2).*T_sq_vector/2);
    tau_sq = 1/gamrnd((p*(p-1)/2+1)/2, 1/rate);    % inv gamma w/ shape=(p*(p-1)/2+1)/2, rate=rate
    xi = 1/gamrnd(1,1/(1+1/tau_sq));    % inv gamma w/ shape=1, rate=1+1/tau_sq

%%% save Omega, lambda_sq, tau_sq
    if iter >burnin           
         omega_save(:,:,iter) = Omega;
		 %omega_vector_save(:,iter) = omega_vector;
         %T_sq_save(:,iter) = T_sq_vector;
         tau_sq_save(1, iter) = tau_sq;
    end

end

end

