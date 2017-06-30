function [results_file,base_error_str] = run_estimators(file_name,alpha)
% function [results_file] = run_estimators(file_name)
%
%Function for running various estimators (S_hat, thresh(S-hat), TME, thresh(TME))
%on general elliptical data contained in ../Data/file_name
%
%Performance is compared against population covariance
%
%Output is [file where data is stored,base_name of files]


disp(file_name)

M = str2num(file_name(2));

load(file_name);

%name of results_file and open (near) empty .mat file
results_file = strcat('Estimators_',file_name);


save(results_file,'ns','gam','iter');

gammaprint=sprintf('%2.1f',gam);
gammaprint = strrep(gammaprint,'.','_');

%matrices containing estimator errors
base_error_str = strcat('gam',gammaprint,'M',num2str(file_name(2)));

E1_errors_str = strcat(base_error_str,'E1');
eval(strcat(E1_errors_str,' = nan(length(ns),iter);'));

E2_errors_str = strcat(base_error_str,'E2');
eval(strcat(E2_errors_str,' = nan(length(ns),iter);'));

E3_errors_str = strcat(base_error_str,'E3');
eval(strcat(E3_errors_str,' = nan(length(ns),iter);'));

E4_errors_str = strcat(base_error_str,'E4');
eval(strcat(E4_errors_str,' = nan(length(ns),iter);'));



for k = 1:length(ns)
    n = ns(k);
    p = round(gam*n);
    nprint=sprintf('%4.4d',n);
    pprint=sprintf('%4.4d',p);
    
    for i = 1:iter
        
        iprint=sprintf('%4.4d',i);
        
        data_str = strcat('M',num2str(M),'p',pprint,'n',nprint,'i',iprint);
        pop_cov_str = strcat('E0',data_str(3:end));
        
        eval(strcat('data = ',data_str,';'));
        eval(strcat('pop_cov = ',pop_cov_str,';'));
        
        %%run various estimators on the data
        
        %Sample Covariance
        scov = data*data'/n;
        
      
         scov = p*scov/trace(scov);
         
         
        
        eval(strcat(E1_errors_str,'(k,i) = norm(scov - pop_cov)/norm(pop_cov);'));
        
        
%         scov_str = strcat('M1',data_str(3:end));
%         eval(strcat(scov_str,' = scov;'));
        
        %Thresholded Sample Covariance
        
        t = sqrt(log(p)/n);
        thresh_scov = hard_thresh(scov,t);
        
        eval(strcat(E2_errors_str,'(k,i) = norm(thresh_scov - pop_cov)/norm(pop_cov);'));
        
        
%         thresh_scov_str = strcat('M2',data_str(3:end));
%         eval(strcat(thresh_scov_str,' = thresh_scov;'));
        
        %TME

        TME = regTME(data,alpha);
        TME = p*TME/trace(TME);
        
        eval(strcat(E3_errors_str,'(k,i) = norm(TME - pop_cov)/norm(pop_cov);'));
        
%         TME_str = strcat('M3',data_str(3:end));
%         eval(strcat(TME_str,' = TME;'));
        
        %thresholded TME
        
        thresh_TME = hard_thresh(TME,t);
        
        eval(strcat(E4_errors_str,'(k,i) = norm(thresh_TME - pop_cov)/norm(pop_cov);'));
        
%         thresh_TME_str = strcat('M4',data_str(3:end));
%         eval(strcat(thresh_TME_str,' = thresh_TME;'));
        
%         %Save files into results matrix
%         save(results_file,pop_cov_str,scov_str,thresh_scov_str,TME_str,thresh_TME_str,'-append')
        
    end
   
    
end

eval(strcat('save(results_file, E1_errors_str,E2_errors_str,E3_errors_str,E4_errors_str,''-append'')'));

