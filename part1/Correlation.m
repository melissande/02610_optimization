function K = Correlation(x)
%Useful constants
m = length(x);

%Create the sum defined by rho and the Threshold
%Rho
rho = ones(1,m-1);
for i = 1:m-1
    rho(i) = x(i)*x(i+1);
end
rho = norm(sum(rho));

%Threshold
T = (1 / sqrt(m-1))*sum(x.^2);

%Test whether the correlation criteria holds
    
if rho > T 
    fprintf('The residuals appear correlated since the autocorrelation p = %f is larger than the threshold T =  %f', rho,T)
    K = true;
else
    fprintf('The residuals appear UNcorrelated since the autocorrelation p = %f is smaller than the threshold T = %f', rho,T)
    K = true;
end
end 