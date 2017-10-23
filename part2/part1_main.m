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
plot(t,y,'o',t,y-r_star,'+')
xlabel('Time in hours')
ylabel('Concentration of NO in ?g/m3 ')
legend('Obersvations','Predictions')

%% Q3: Optimal order of the fit

%% Test for random sign
n=3;
z_plus_moins=2;
while z_plus_moins>=1.96
    n=n+2;
    if n>24
        
        disp ('stop')
        fprintf('z+-=%d',z_plus_moins)
        break
    end
    [x_star,r_star] = NOfit(t,y,n);
    res=r_star;
    res(res==0)=[];
    m=length(res);
    bin_sign_r=res>0;
    N_runs=nnz(diff(bin_sign_r))+1;
    n_plus=sum(bin_sign_r);
    n_moins=sum(~bin_sign_r);
    mu_u=(2*n_plus+n_moins)/m+1;
    var_u=(mu_u-1)*(mu_u-2)/(m-1);
    z_plus_moins=abs(N_runs-mu_u)/sqrt(var_u);
    
end

norm_r_star=sqrt(sum((r_star).^2));
fprintf('The norm of the residual is %d  and the order of the model is %d\n',norm_r_star,n)

figure(1);
plot(t,y,'o',t,y-r_star,'+')
xlabel('Time in hours')
ylabel('Concentration of NO in microg/m3 ')
legend('Obersvations','Predictions')
%% Test for autocorrelation

n=3;
rho=3;
t_rho=2;

while abs(rho)>t_rho
    
    [x_star,r_star] = NOfit(t,y,n);
    res=r_star;

    ri=res;
    ri(end)=[];
    riplus1=res;
    riplus1(1)=[];
    rho=sum(ri.*riplus1);
    t_rho=1/sqrt(size(res,1)-1)*sum(res.^2);
    n=n+2;
end

norm_r_star=sqrt(sum((r_star).^2));
fprintf('The norm of the residual is %d  and the order of the model is %d\n',norm_r_star,n)

figure(1);
plot(t,y,'o',t,y-r_star,'+')
xlabel('Time in hours')
ylabel('Concentration of NO in microg/m3 ')
legend('Obersvations','Predictions')

%% Estimating the Standard Deviation of the Solution Coefficients
m=length(t);
n=zeros((m/2-1),1);
s=zeros((m/2-1),1);
n(1)=3;
k=1;
while n(k)<=m
    [x_star,r_star] = NOfit(t,y,n(k));
    s(k)=sqrt(sum(r_star.*r_star)/(m-n(k)));
    n(k+1)=n(k)+2;
    k=k+1;
end

figure(2)
plot(n(1:end-1),s,'o')
xlabel('Order n')
ylabel('Scaled residual norm s')

    
%n=19 ->solution