%% PROBLEM 1

%% Q1: PLOT DATA

data_pb1=load('Problem1Data.mat');

ti=data_pb1.t;
y_obs=data_pb1.y;

alpha=1.0;
beta=0.0;

y_true=alpha*ti+beta;

figure(1)
plot(ti,y_obs,'bo',ti,y_true,'r');
legend('Observations','True Model')
xlabel('t')
ylabel('y')

%% Used in all problems
n=length(y_true);
A=[ti,ones(n,1)];
b=y_obs;
p=2;


%% Q2: l2 Squares Estimation

H=A'*A;
g=-A'*b;
gamma=-1/2*b'*b;

x_ls=-inv(H)*g;

y_ls=A*x_ls;
figure(2)
plot(ti,y_obs,'bo',ti,y_true,'r',ti,y_ls,'g');
title('Least Square Solution')
legend('Observations','True Model','Least Square Solution')
xlabel('t')
ylabel('y')

%Histogram with outliers
res_ls=y_obs-y_ls;
figure(3)
subplot(1,2,1)
hist(res_ls,50);
title('Histogram of the errors for Least Square Solution')
%Histogram without outliers
idx=find(res_ls>12);
res_ls_wo=res_ls;
res_ls_wo(idx)=[];
n2=length(res_ls_wo);



subplot(1,2,2)
hist(res_ls_wo,50);
title('Histogram of the errors for Least Square Solution without outliers')

%Goodness of the fit

%Predictions
std_noise=sqrt(sum(res_ls_wo.^2)/length(res_ls_wo));
PI_data_ls=y_ls+tinv([0.025  0.975],n2-p)*std_noise;
figure(4)
plot(ti,y_obs,'bo',ti,y_true,'r',ti,y_ls,'g',ti,PI_data_ls(:,1),'--g',ti,PI_data_ls(:,2),'--g');
legend('Observations','True Model','Least Square Solution','Low Bound PI',' High Bound PI')
xlabel('t')
ylabel('y')

%Parameters
cov_param=cov(x_ls);
std_param=diag(sqrt(cov_param));
PI_param_ls=x_ls+tinv([0.025  0.975],p)*std_param;
disp('LS: PI for parameters')
disp(PI_param_ls)


%% Q3: l1  Estimation

f_lp=[zeros(p,1);ones(n,1)];
A_lp=[-A,-eye(n);A,-eye(n)];
b_lp=[-b;b];
  
x_l1 = linprog(f_lp,A_lp,b_lp); 
x_l1=x_l1(1:2);
y_l1=A*x_l1;

figure(5)
plot(ti,y_obs,'bo',ti,y_true,'r',ti,y_l1,'g');
legend('Observations','True Model','L1 Estimate')
title('L1 Estimate')
xlabel('t')
ylabel('y')

%Histogram with outliers
res_l1=y_obs-y_l1;
figure(6)
subplot(1,2,1)
hist(res_l1,50);
title('Histogram of the errors for L1 Estimate Solution')
%Histogram without outliers
idx=find(res_l1>12);
res_l1_wo=res_l1;
res_l1_wo(idx)=[];
n2=length(res_l1_wo);


subplot(1,2,2)
hist(res_l1_wo,50);
title('Histogram of the errors for L1 Estimate without outliers')

%Goodness of the fit

%Predictions
std_noise=sqrt(sum(res_l1_wo.^2)/length(res_l1_wo));
PI_data_l1=y_l1+tinv([0.025  0.975],n2-p)*std_noise;

figure(7)
plot(ti,y_obs,'bo',ti,y_true,'r',ti,y_l1,'g',ti,PI_data_l1(:,1),'--g',ti,PI_data_l1(:,2),'--g');
legend('Observations','True Model','L1 Estimate','Low Bound PI',' High Bound PI')
xlabel('t')
ylabel('y')

%Parameters
cov_param=cov(x_l1);
std_param=diag(sqrt(cov_param));
PI_param_l1=x_l1+tinv([0.025  0.975],p)*std_param;
disp('L1: PI for parameters')
disp(PI_param_l1)

