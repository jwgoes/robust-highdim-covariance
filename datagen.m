function [data_file] = datagen(ns,gam,iter,rho,u)
% function [data_file] = datagen(ns,gam,iter,rho,u)
%
% Script for generating data from a generalized elliptical model. Output is
% file where data is stored.
% 
% Data is drawn from a sparse covariance matrix with parameter rho,
% described in Bickel and Levina's paper
% 
% Generate (p,n) pairs
% where (p,n) is determined by a vector of ns
% and the paramter gam = p/n
%
%iter = number of iterations
% 
% 
% %scale Gaussian data by u
% u = 'cauchy'; M=1;
% %u = 'laplace'; M=2;
% %u = 'normal'; M=3;

%stringify gam
gammaprint=sprintf('%2.1f',gam);
 gammaprint = strrep(gammaprint,'.','_');

if strcmp(u,'cauchy') == 1
    M=1;
elseif strcmp(u,'laplace') == 1
    M=2;
elseif strcmp(u,'normal') == 1
    M=3;
else
    error('not a valid u option')
end


%Name output .mat file of results
data_file = strcat('M',num2str(M),'gam',gammaprint,'n',num2str(min(ns)),'to',num2str(max(ns)));


save(data_file,'ns','gam','iter');




for i = 1:iter
    
    for n = ns
        
        p = round(gam*n); %p is determined by n and gamma
        
        %the following puts n,p and i in the form of a 4-digit string
        nprint=sprintf('%4.4d',n);
        pprint=sprintf('%4.4d',p);
        iprint=sprintf('%4.4d',i);
        
        % create Bickel and Levina's sparse covariance matrix
        exponents = [0:p-1];
        one_row = rho.^exponents;
        covv = toeplitz(one_row);
        
        %name population covariance
        truth_name = sprintf('E0p%sn%si%s',pprint,nprint,iprint);
        
        
        eval(strcat(truth_name,' = covv;'));
        
        
        %generate Gaussian data, N(0,covv)
        
        data = sqrtm(covv)*randn(p,n);
        

        
        
        % scale data by u_vec
        if strcmp(u,'cauchy') == 1
            u_vec = tan((rand(n,1)-.5)*pi);
            M=1;
        elseif strcmp(u,'laplace') == 1
            unif = rand(n,1)-.5;
            u_vec = sign(unif).*log(1-2*abs(unif));
            M=2;
        elseif strcmp(u,'normal') == 1
            u_vec = ones(n,1);
            M=3;
        else
            error('Not a valid option for u-scaling')
        end
        
        %Elliptical data
        data = data*diag(u_vec);
        
        
        %name data matrix
        data_name = sprintf('M%sp%sn%si%s',num2str(M),pprint,nprint,iprint);
        
        eval(strcat(data_name,' = data;'));
        
        save(data_file,truth_name,data_name,'-append'); 
        
    end
    
    
    
end



