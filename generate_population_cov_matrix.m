function Sp = generate_population_cov_matrix(p,rho)
%function Sp = generate_population_cov_matrix(p,rho)
% create Bickel and Levina's sparse covariance matrix


Sp = zeros(p,p); 
  
exponents = [0:p-1];
one_row = rho.^exponents;
Sp = toeplitz(one_row);
  