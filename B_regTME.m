function [cov,weights,iter]=B_regTME(X,alpha,epsilon,max_iter)
%function [cov,weights,iter]=B_regTME(X,alpha,epsilon,max_iter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularized Tyler M-Estimator (TME)
%
% X: data matrix with each column representing a point
% alpha: regularization coefficient (default provided)
% cov: the estimated covariance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %fprintf('inside B_regTME nargin %d\n',nargin); 
 
if nargin < 4 
    max_iter= 2000;
end
if nargin < 3
    epsilon = 10^(-5); 
end
if nargin==1
    alpha=min(1,p/n);  %default alpha. Existence results require alpha>p/n-1
end


format long

[p,n]=size(X);

initcov=eye(p); %Initialize with identity matrix - 
                %any positive definite initializations is possible

%allocate variables used in iterations
oldcov=initcov+1; 
cov=initcov;
iter=0;
weights=zeros(n,1);


while (norm((p*oldcov/trace(oldcov)-p*cov/trace(cov)),'fro')>epsilon & iter<max_iter)

    oldcov=cov;
    cov=cov+alpha*eye(p)/(1+alpha);
    
    weights=sum(X.*(cov\X))'; 
    weights=1./weights;  
    
    cov=(X.*repmat(((weights)),1,p)')*X';
    cov=.5*(cov+cov');   % forces matrix to be symmetric
    
    cov = (p/n)*cov/(1+alpha);
    iter=iter+1;

end


%fprintf('exit B_regTME iter %d\n',iter); 




