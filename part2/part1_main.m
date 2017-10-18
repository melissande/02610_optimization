%% PART 1

t=(1:24)';
y=[91.22,28.04,22.91,26.65,42.96,101.05,202.36,328.02,364.12,299.23,238.00,...
    227.49,218.03,223.62,238.75,271.26,267.72,251.32,230.04,...
    206.69,170.77,131.76,143.85,157.57]';

%% Q1: Matlab Code
%See NOfit.m
%% Q2: Test your software
n=3;
[x_star,r_star] = NOfit(t,y,n);

norm_r_star=sqrt(sum((r_star).^2));
fprintf('The model vector is [%d,%d,%d] and the norm of the residual is %d \n',x_star,norm_r_star)

figure(1);
plot(t,y,'o',t,y-r_star,'-')
xlabel('Time in hours')
ylabel('Concentration of NO in ?g/m3 ')
legend('Obersvations','Predictions')

%% Q3: Optimal order of the fit
m=length(r_star);
sign_r=sign(r_star);
n_plus=sum(sign_r>0);
n_moins=sum(sign_r<0);
mu_u=(2*n_plus+n_moins)/m+1;
var_u=(mu_u-1)*(mu_u-2)/(m-1);