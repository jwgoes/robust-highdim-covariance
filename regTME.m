function [cov]=regTME(X,alpha)
%function [cov]=regTME(X,alpha)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularized Tyler M-Estimator (TME)
%
% X: data matrix with each column representing a point
% alpha: regularization coefficient (default provided)
% cov: the estimated covariance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long

[p,n]=size(X);
if nargin==1
    alpha=min(1,p/n);  %default alpha. Existence results require alpha>p/n-1
end

cov=eye(p); %Initialize with identity matrix - 
                %any positive definite initializations is possible

%allocate variables used in iterations
iter=0;
weights=zeros(n,1);


err=10;

while (err>10^-12 & iter<1400)

    oldcov=cov;
    cov=cov+alpha*eye(p)/(1+alpha);
    
    weights=sum(X.*(cov\X))'; 
    weights=1./weights;  
    
    cov=(X.*repmat(((weights)),1,p)')*X';
    cov=.5*(cov+cov');
    
    cov = (p/n)*cov/(1+alpha);
    iter=iter+1;
    
    err = norm((p*oldcov/trace(oldcov)-p*cov/trace(cov)),'fro');

end



