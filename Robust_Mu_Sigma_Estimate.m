function [mu,sigma] = Robust_Mu_Sigma_Estimate(weights,r)
%function [mu,sigma] = Robust_Mu_Sigma_Estimate(weights,r)

weights_s = sort(weights); 
n = length(weights); 

xi = linspace(weights_s(1),weights_s(end),500); 

[f] = ksdensity(weights_s,xi); 

[f_max idx] = max(f); 

mu = xi(idx); 

% finding the right edge where f > r * f_max
for j=idx+1:n
    if f(j) < r* f_max
        break; 
    end
end
idxR = j-1; 
xR = xi(idxR);

for j=idx-1:-1:1
    if f(j)< r* f_max
        break; 
    end
end
idxL = j+1; 
xL = xi(idxL); 

x_minus_mu = 1/2 * abs(xR-mu) + 1/2 * abs(xL-mu); 
sigma = x_minus_mu*sqrt(-0.5/log( r ) ) ; 

if 0 
	figure(40); clf;  plot(xi,f,'b.-'); grid on; title(['KERNEL DENSITY ESTIMATE \mu = ' num2str(mu) ' \sigma = ' num2str(sigma) ]); 
	hold on; 
	plot(xi(idxL:idxR),f(idxL:idxR),'go'); 

	figure(11); clf; plot(1:n,weights,'ko-','markersize',8); grid on; 
	drawnow; 
end
%pause; 

return; 

