% Measure sensitivity of runtime and performance of regularized-TME
% on the parameter alpha
%  The goal of this script is to observe how the choice of regularization
%  constant is affected by a changing ratio of p/n
%
%  The dimension 'p' is fixed at a chosen value. 'ns is the set of n, the
%  number of samples. 'alphas' determines the EXACTLY 4 regularization
%  constants to try. 'iter' is the number of iterations to run.
%
%  Results are automatically plotted at the end of the experiment
%%%%%%%%%%%%%%%
close all
clear all

iter = 100;

pns = [[800 400]' [800 200]' [ 400 200]']; 





alphas = [[.1:.1:.9],[1:20]];
alphas = [[.2:.2:.8],[1:20]];





%loop over and (p,n), alphas, iter

alpha_times = nan(iter,length(alphas),3);
alpha_perfs = nan(iter,length(alphas),3);

SpNorms = nan(3,1);

for m=[1:3]


p = pns(1,m);
n = pns(2,m);


%generate BL covariance
exponents = [0:p-1];
one_row = .7.^exponents;
Sp = toeplitz(one_row);

SpNorms(m) = norm(Sp);


for i = [1:iter]

iter_p_n = [i,p,n]
    
    %generate data
    
    data = sqrtm(Sp)*randn(p,n);
    
    for k=[1:length(alphas)]
        alpha = alphas(k);
        
        [i,k/length(alphas)]
        
        %run TME while recording runtime
        timer1 = tic;
        TME = regTME(data,alpha);
        
        alpha_times(i,k,m) = toc(timer1);
        TME = p*TME/trace(TME);
        TME = hard_thresh(TME,sqrt(log(p)/n));
        
        if trace(TME) == 0
            error('problem')
        end


		err = norm(TME - Sp)/norm(Sp);

        alpha_perfs(i,k,m) = err;
        
        
        %         alpha_perfs(i,k) = norm(p*TME/trace(TME) - Sp);
        
    end
    
    
end

end

%save('alpha_data3.mat','alpha_times','alpha_perfs')


% %plot
%
alpha_times_mean = log(mean(alpha_times));
alpha_perfs_mean = log(mean(alpha_perfs));













apm1 = reshape(alpha_perfs_mean(1,:,1),24,1);
atm1 = reshape(alpha_times_mean(1,:,2),24,1);

apm2 = reshape(alpha_perfs_mean(1,:,1),24,1);
atm2 = reshape(alpha_times_mean(1,:,2),24,1);


apm3 = reshape(alpha_perfs_mean(1,:,1),24,1);
atm3 = reshape(alpha_times_mean(1,:,2),24,1);




figure;









sz=200;
c=[1 .5 0];

ppp=plot(alphas,apm1/SpNorms(1),'^',alphas,apm2/SpNorms(2),'o',alphas,apm3/SpNorms(3),'s');
set(ppp(1),'MarkerEdgeColor','None','MarkerFaceColor','blue','MarkerSize',24);
set(ppp(2),'MarkerEdgeColor','None','MarkerFaceColor',c,'MarkerSize',24);
set(ppp(3),'MarkerEdgeColor','None','MarkerFaceColor','black','MarkerSize',24);
xlabel('\alpha');
ylabel('LRE');
set(gca, 'FontSize', 36)
h = legend('p=800,n=200','p=400,n=200','p=800,n=400');
set(h,'FontSize',36);

figure;
ppp = plot(alphas,atm1/SpNorms(1),'^',alphas,atm2/SpNorms(2),'o',alphas,atm3/SpNorms(3),'s')
set(ppp(1),'MarkerEdgeColor','None','MarkerFaceColor','blue','MarkerSize',24);
set(ppp(2),'MarkerEdgeColor','None','MarkerFaceColor',c,'MarkerSize',24);
set(ppp(3),'MarkerEdgeColor','None','MarkerFaceColor','black','MarkerSize',24);
xlabel('\alpha');
ylabel('log(seconds)');
set(gca, 'FontSize', 36)
h = legend('p=800,n=200','p=400,n=200','p=800,n=400');
set(h,'FontSize',36);

