%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script that produces creates data, runs estimators on the data and plots
% results
%
% % datagen(ns,gam,iter,rho,u) produces data with parameters:
% % ns - vector of ns
% % gam = gamma = p/n
% % iter - number of iterations
% % rho - parameter rho controlling sparsity of BL-covariance
% % u in {'cauchy','laplace','normal'}
%
% % run_estimators(-) runs the various estimators on the data produced
% above
%
% % the data is then plotted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%parameters for data generation
ns = [200:50:400];

gamma = .5; %gamma = p/n determines a p for every n in ns

iter = 20; %number of iterations 

rho = .7; %

%hu = 'cauchy';
u = 'laplace';
 %u = 'normal';

alpha=10; %regularization parameter

%%data generation
[data_file] = datagen(ns,gamma,iter,rho,u);


%%run estimators

[results_file,base_name] = run_estimators(data_file,alpha);


%%load results



load(results_file)

eval(strcat('E1 = ',base_name,'E1'));
eval(strcat('E2 = ',base_name,'E2'));
eval(strcat('E3 = ',base_name,'E3'));
eval(strcat('E4 = ',base_name,'E4'));


size(E1)

EE1 = log(mean(E1,2));
EE2 = log(mean(E2,2));
EE3 = log(mean(E3,2));
EE4 = log(mean(E4,2));





left_color = [0 0 0];
right_color = [0 0 0];
figure;
ppp=plot(ns,EE1',':bs', ns,EE2','--mo', ns,EE3',':rs', ns,EE4','-ko');
legend('SampCov','th-SampCov','TME','th-TME');
yyaxis left;
xlabel('n');
ylabel('LRE');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(ppp(1),'MarkerEdgeColor','None','MarkerFaceColor','b','MarkerSize',12)
set(ppp(2),'MarkerEdgeColor','b','MarkerFaceColor','None','MarkerSize',12)
set(ppp(3),'MarkerEdgeColor','None','MarkerFaceColor','r','MarkerSize',12)
set(ppp(4),'MarkerEdgeColor','r','MarkerFaceColor','None','MarkerSize',12)
set(gca, 'FontSize', 16)
%xlim([100,1000]);
title(strcat('u_i\sim ',u,'(0,1)'))
ax = gca;
ax.YColor = left_color;
yyaxis right
axis tight
%set(gca,'xtick',[],'ytick',[],'layer','bottom','box','on')
ylabel('\gamma = .5','rot',270,'VerticalAlignment','bottom','FontSize', 24)
ax = gca;
ax.YColor = right_color;






