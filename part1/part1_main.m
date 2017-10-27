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

m=length(t);
%Used to check the algo is working but it's not .. for the moment!
res=zeros(m,(m/2-1)); %residuals for each order
Bin=zeros(m,(m/2-1)); %Binary (1 if res>0 and 0 if res<0)
N_runs=zeros((m/2-1),1);%nb of run for each ite
N_plus=zeros((m/2-1),1); % nb of + for each iteration
N_moins=zeros((m/2-1),1);% nb of - for each iteration
%%%%
z=zeros((m/2-1),1);
order=zeros((m/2-1),1);
n=3;

z_plus_moins=2;
k=1;
while n<=24

    
    [x_star,r_star] = NOfit(t,y,n);
    res(:,k)=r_star;
    %res(res==0)=[];

    bin_sign_r=res(:,k)>0;
    Bin(:,k)=bin_sign_r;
    n_runs=nnz(diff(bin_sign_r))+1;
    n_plus=sum(bin_sign_r);
    n_moins=sum(~bin_sign_r);
    mu_u=(2*n_plus+n_moins)/m+1;
    var_u=(mu_u-1)*(mu_u-2)/(m-1);
    z_plus_moins=abs(n_runs-mu_u)/sqrt(var_u);
    
    N_runs(k)=n_runs;
    N_plus(k)=n_plus;
    N_moins(k)=n_moins;
    z(k)=z_plus_moins;
    order(k)=n;
    k=k+1;
    n=n+2;
end

%TO FIND orders that work, Have a look to order() and  corresponding z()
%% Test for autocorrelation
m=length(t);
ratio=zeros((m/2-1),1);
order=zeros((m/2-1),1);

k=1;
n=3;
rho=3;
t_rho=2;

while n<=24
    
    [x_star,r_star] = NOfit(t,y,n);
    res=r_star;

    ri=res;
    ri(end)=[];
    riplus1=res;
    riplus1(1)=[];
    rho=sum(ri.*riplus1);
    t_rho=1/sqrt(size(res,1)-1)*sum(res.^2);
    order(k)=n;
    ratio(k)=abs(rho)/t_rho;
    n=n+2;
    k=k+1;
end

%n=7-> min order solution
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

%n=7 ->WORKS!
n=7; %which is s(3)
[x_star,r_star,A] = NOfit(t,y,n);
sigma=s(3);
cov_x=sigma^2*inv(A'*A);
std_x=sqrt(diag(cov_x));

