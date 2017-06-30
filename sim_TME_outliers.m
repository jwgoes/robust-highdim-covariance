%n = 400; p = 800; 

n = 500; p = 500; 
alpha = 4;

n_sigma = 2.0; 

rho = 0.7; 

iter = 10; 
seed = 1; randn('state',seed); 

outlier_ratio = [0:0.05:0.5]; 
%outlier_ratio = [0.4:0.05:0.5]; 
n_epsilon = length(outlier_ratio); 

Sp = generate_population_cov_matrix(p,rho); 
Sp = p*Sp/trace(Sp); 
norm_Sp = norm(Sp); 
Sp_sqrt = sqrtm(Sp); 

% generate outlier covariance
if 1
    Sout = diag(4*rand(p,1)+1);
else
    Sout = eye(p);
    %Sout(1,1) = sqrt(p);  Sout(2,2) = sqrt(p); 
    Sout(1,1) = p; Sout(2,2) = p/2; 
end

Sout = p * Sout/trace(Sout); 
Sout_sqrt = sqrtm(Sout); 

Xtmp = rand(p,p); [U,sigma,V] = svd(Xtmp); 

%%
err_Sn = zeros(iter,n_epsilon); 
err_TME= zeros(iter,n_epsilon);
err_TME_th=zeros(iter,n_epsilon); 
err_Sn_th = zeros(iter,n_epsilon);
err_Out_TME=zeros(iter,n_epsilon); 
err_Out_TME_th=zeros(iter,n_epsilon); 

W = zeros(n_epsilon,n); 

mu_sigma = zeros(n_epsilon,2); 

for counter = 1:iter
    counter
    Xin = Sp_sqrt * randn(p,n);     %X of size p x n 
    
    Xout = V * Sout_sqrt * V' * randn(p,n);   %Xout also of size p x n 
    
    norm_arr = zeros(2*n,1); 
    for j=1:n
        norm_arr(j) = norm(Xin(:,j)); 
        norm_arr(n+j)=norm(Xout(:,j)); 
    end
    %figure(44); plot(norm_arr/sqrt(p),'b.-');    return; 
    for j=1:n
        Xin(:,j) = Xin(:,j)/norm(Xin(:,j));
        Xout(:,j) = Xout(:,j)/norm(Xout(:,j)); 
    end
    
    for j=1:n_epsilon
        
        n_out = floor(outlier_ratio(j) * n);
        
        X = Xin; 
        X(:,1:n_out) = Xout(:,1:n_out); 
        
        Sn = 1/n * X * X';         % sample covariance matrix 
        Sn = 1/2 * (Sn + Sn'); 
    
        Sn = p*Sn/trace(Sn); 
        
         
        threshold_TME = 1e-12; 
        max_iter_TME=1000; 
        [TME_hat,TME_w,tme_iter] = B_regTME(X,alpha,threshold_TME,max_iter_TME); 
    
        
        W(j,:) = TME_w'; 
        r = 0.7; 
        [mu_weight,sigma_weight]= Robust_Mu_Sigma_Estimate(TME_w,r); 
        
        mu_sigma(j,:) = [mu_weight, sigma_weight]; 
        
        figure(2); clf;
        plot(1:n,TME_w,'kp--'); grid on; hold on; 
        plot(1:n,mu_weight*ones(n,1),'b-',1:n,(mu_weight-n_sigma*sigma_weight)*ones(n,1),'g.-','linewidth',2);
        plot(1:n,(mu_weight+n_sigma*sigma_weight)*ones(n,1),'g.-','linewidth',2); 
        title(['TME weights epsilon ' num2str(outlier_ratio(j)) ]); 
        drawnow; 
        
        %% INLIERS
        th_weight = TME_w .*( abs(TME_w - mu_weight)<n_sigma*sigma_weight); 
        idx_in = find(th_weight>0); 
        n_in = sum(th_weight>0);
        
        [Out_TME, new_weights] = B_regTME(X(:,idx_in),alpha,threshold_TME,max_iter_TME);
        %Out_TME = X.*(repmat(th_weight,1,p)')*X'; 
        
        figure(55); plot(new_weights,'rs-'); grid on; title('NEW WEIGHTS\n'); drawnow; 
        
        Out_TME = p * Out_TME / trace(Out_TME); 
        
 
        fprintf('n_in %d\n',n_in); 
        
        Normalized_TME_hat = p * TME_hat / trace(TME_hat); 
    
        err_Sn(counter,j) = norm(Sn - Sp); 
        err_TME(counter,j)= norm(Normalized_TME_hat-Sp);
        err_Out_TME(counter,j) = norm(Out_TME-Sp); 
        
        threshold_out = sqrt(log(p)/n_in); 
        Out_TME_th = Out_TME .* (abs(Out_TME) > threshold_out) ; 
        err_Out_TME_th(counter,j) = norm(Out_TME_th-Sp); 
        
        threshold = sqrt(log(p)/n); 
        Sn_th = Sn .* (abs(Sn) > threshold);
    
        TME_hat = Normalized_TME_hat .* ( abs(Normalized_TME_hat) > threshold); 
    
        err_Sn_th(counter,j) = norm(Sn_th-Sp); 
        err_TME_th(counter,j)= norm(TME_hat-Sp); 
    end
    if mod(counter,2)==0
 %%       
        figure(42); clf; set(gca,'fontsize',20); 
        plot(outlier_ratio,log(mean(err_Sn(1:counter,:))/norm_Sp),'b.-',outlier_ratio,log(mean(err_TME(1:counter,:))/norm_Sp),'rs-'); 
        hold on; 
        plot(outlier_ratio,log(mean(err_Sn_th(1:counter,:))/norm_Sp),'k>--',outlier_ratio,log(mean(err_TME_th(1:counter,:))/norm_Sp),'gd-'); 
        plot(outlier_ratio,log(mean(err_Out_TME(1:counter,:))/norm_Sp),'cp-',outlier_ratio,log(mean(err_Out_TME_th(1:counter,:))/norm_Sp),'md-'); 
        title(['n = ' num2str(n) ' p= ' num2str(p)  ' iter= ' num2str(counter)]); 
        xlabel('outlier ratio'); 
        legend('Scov','TME','th(Scov)','th(TME)','TME(outliers removed)','th(TME w/o outliers)'); 
        grid on ;
        drawnow; 
    end
end


%%
        figure(42); clf; set(gca,'fontsize',28); 
        hold on; 
        
        plot(outlier_ratio,log(mean(err_TME_th(1:counter,:))/norm_Sp),'bo-','markersize',12,'linewidth',2); 
        plot(outlier_ratio,log(mean(err_Out_TME_th(1:counter,:))/norm_Sp),'md-','markersize',12,'linewidth',2); 
        title(['$n$ = ' num2str(n) ' $p$= ', num2str(p)  ', $\alpha$= ' num2str(alpha)],'Interpreter','latex'); 
        xlabel('outlier ratio'); 
        ylabel('LRE'); 
        axis([0 0.45 min(log(mean(err_TME_th(1:counter,:))/norm_Sp))*1.1 max(log(mean(err_Out_TME_th(1:counter,end-1))/norm_Sp),log(mean(err_TME_th(1:counter,end-1))/norm_Sp))*1.1]); 
        legend('th-RegTME','Outliers Removed','Location','northwest'); 
        %legend('Scov','TME','th(Scov)','th(TME)'); 
        grid on ; box on; 
        
        drawnow; 

        
 %%
 j = find(outlier_ratio==0.2); 
 figure(33); clf; set(gca,'fontsize',28); 
        hold on; 
 plot(W(j,:),'ko'); grid on; 
 mu = mu_sigma(j,1); sigma = mu_sigma(j,2); 
 plot(1:n,mu*ones(n,1),'b-',1:n,(mu-n_sigma*sigma)*ones(n,1),'r-.','linewidth',3);
        plot(1:n,(mu+n_sigma*sigma)*ones(n,1),'r-.','linewidth',3);
        xlabel('$j$','interpreter','latex'); ylabel('$w_j$','interpreter','latex'); 
        
        axis([1 n min(W(j,:))*0.95 max(W(j,:))*1.05]); 
        title(['TME weights $\epsilon=$ ' num2str(outlier_ratio(j)) ],'interpreter','latex'); 
        box on; 

 %%
  j = find(outlier_ratio==0.4); 
 figure(55); clf; set(gca,'fontsize',28); 
        hold on; 
 plot(W(j,:),'ko'); grid on; 
 mu = mu_sigma(j,1); sigma = mu_sigma(j,2); 
 plot(1:n,mu*ones(n,1),'b-',1:n,(mu-n_sigma*sigma)*ones(n,1),'r-.','linewidth',3);
        plot(1:n,(mu+n_sigma*sigma)*ones(n,1),'r-.','linewidth',3);
        xlabel('$j$','interpreter','latex'); ylabel('$w_j$','interpreter','latex'); 
        title(['TME weights $\epsilon=$ ' num2str(outlier_ratio(j)) ],'interpreter','latex'); 
        axis([1 n min(W(j,:))*0.95 max(W(j,:))*1.05]);
        
        axis([1 n min(W(j,:))*0.95 7.5]);

        box on; 

        