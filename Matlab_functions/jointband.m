function [upper_CI, lower_CI] = jointband(MCMC_P,alpha)
  % function to calculate (1-alpha) joint credible band 
  
   % Inputs
   %    'MCMC_P' - a B by T matrix containing MCMC samples. 
   %        B - number of MCMC iterations (MCMCspecs.B)
   %        T - number of function samples (number of columns in Y)
   %    'alpha' - a row vector containing significance levels at which 
   %        joint credible bands are to be returned. 
   %    'region' - locations in T where the joint bands are computed
   
   % Outputs
   %   'upper_CI' - a length(alpha) by length(region) matrix containing the upper bounds 
   %        of the joint credible bands. The first row corresponds to 
   %        the first level in alpha.       
   %   'lower_CI' - a length(alpha) by length(region) matrix containing the lower bounds 
   %        of the joint credible bands. 
    
   [B,T] = size(MCMC_P); 
   sd_P = NaN(1,T);
   for i=1:T
       sd_P(i) = std(MCMC_P(:,i));
   end
   mean_P = mean(MCMC_P);
   
   z_P = NaN(1,B);
   for j=1:B
        z_P(j) = max(abs((MCMC_P(j,:)-mean_P)./sd_P));
   end
   
   c = quantile(z_P,1-alpha);
   upper_CI = mean_P + c*sd_P;
   lower_CI = mean_P - c*sd_P;
   
end

