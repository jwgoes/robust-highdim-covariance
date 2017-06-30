%Measure sensitivity of runtime and performance of regularized-TME
%on the parameter alpha
%  The goal of this script is to observe how the choice of regularization
%  constant is affected by a changing ratio of p/n
%
%  The dimension 'p' is fixed at a chosen value. 'ns is the set of n, the
%  number of samples. 'alphas' determines the EXACTLY 4 regularization 
%  constants to try. 'iter' is the number of iterations to run.
%
%  Results are automatically plotted at the end of the experiment
%%%%%%%%%%%%%%%

iter = 100;

p = 120;
ns = [5:5:100]
alphas = [1,2,3,4];

p=120;
iter=5;
ns = [20:20:100];

%loop over and (p,n), alphas, iter

alpha_times = nan(length(alphas),length(ns),iter);
alpha_perfs = nan(length(alphas),length(ns),iter);

%generate BL covariance
exponents = [0:p-1];
one_row = .7.^exponents;
Sp = toeplitz(one_row);

for k = 1:length(alphas)
    alpha = alphas(k);
    
    for j = 1:length(ns)
        n = ns(j);
        
        
        
        for i = [1:iter]
            
            iter_nIdx_AlphaIdx = [i,j,k]
            
            %generate data
            
            data = sqrtm(Sp)*randn(p,n);
            
            %run TME while recording runtime
            timer1 = tic;
            TME = regTME(data,alpha);
            alpha_times(k,j,i) = toc(timer1);
            
            alpha_perfs(k,j,i) = norm(p*TME/trace(TME) - Sp)/norm(Sp);
            
        end
        
    end
    
end

%save('alpha_data.mat','alpha_times','alpha_perfs')


% %plot
alpha1_times = log(mean(reshape(alpha_times(1,:,:),length(ns),iter),2));
alpha1_perfs = log(mean(reshape(alpha_perfs(1,:,:),length(ns),iter),2));

alpha2_times = log(mean(reshape(alpha_times(2,:,:),length(ns),iter),2));
alpha2_perfs = log(mean(reshape(alpha_perfs(2,:,:),length(ns),iter),2));

alpha3_times = log(mean(reshape(alpha_times(3,:,:),length(ns),iter),2));
alpha3_perfs = log(mean(reshape(alpha_perfs(3,:,:),length(ns),iter),2));

alpha4_times = log(mean(reshape(alpha_times(4,:,:),length(ns),iter),2));
alpha4_perfs = log(mean(reshape(alpha_perfs(4,:,:),length(ns),iter),2));

%save('test2.mat','alpha_times','alpha_perfs')





figure;
ppp = plot(ns,alpha1_perfs/norm(Sp),'bs',ns,alpha2_perfs/norm(Sp),'mo',ns,alpha3_perfs/norm(Sp),'rs',ns,alpha4_perfs/norm(Sp),'ko')
xlabel('n');
ylabel('LRE');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(ppp(1),'MarkerEdgeColor','None','MarkerFaceColor','b','MarkerSize',36)
set(ppp(2),'MarkerEdgeColor','b','MarkerFaceColor','None','MarkerSize',36)
set(ppp(3),'MarkerEdgeColor','None','MarkerFaceColor','r','MarkerSize',36)
set(ppp(4),'MarkerEdgeColor','r','MarkerFaceColor','None','MarkerSize',36)
xlim([min(ns),max(ns)]);
h_legend = legend(strcat('\alpha = ',num2str(alphas(1))),strcat('\alpha = ',num2str(alphas(2))),strcat('\alpha = ',num2str(alphas(3))),strcat('\alpha = ',num2str(alphas(4))));
set(gca, 'FontSize', 36)
set(h_legend,'FontSize',36);


figure;
ppp=plot(ns,alpha1_times,'bs',ns,alpha2_times,'mo',ns,alpha3_times,'rs',ns,alpha4_times,'ko')
xlabel('n');
ylabel('log(seconds)');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(ppp(1),'MarkerEdgeColor','None','MarkerFaceColor','b','MarkerSize',36)
set(ppp(2),'MarkerEdgeColor','b','MarkerFaceColor','None','MarkerSize',36)
set(ppp(3),'MarkerEdgeColor','None','MarkerFaceColor','r','MarkerSize',36)
set(ppp(4),'MarkerEdgeColor','r','MarkerFaceColor','None','MarkerSize',36)
set(gca, 'FontSize', 36)
xlim([min(ns),max(ns)]);
%lgd2.Location = 'southeast';
h_legend = legend(strcat('\alpha = ',num2str(alphas(1))),strcat('\alpha = ',num2str(alphas(2))),strcat('\alpha = ',num2str(alphas(3))),strcat('\alpha = ',num2str(alphas(4))));
set(h_legend,'FontSize',36);