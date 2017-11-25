%% PROBLEM 1
close all
clear all
clc

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
gamma=-1/2*(b)'*b;

x_ls=-H\g;

y_ls=A*x_ls;
figure(2)
plot(ti,y_obs,'bo',ti,y_true,'r',ti,y_ls,'g');
title('Least Square Solution with Outliers')
legend('Observations','True Model','Least Square Solution')
xlabel('t')
ylabel('y')

%Histogram with outliers
res_ls=y_obs-y_ls;
figure(3)
subplot(1,2,1)
hist(res_ls,50);
title({'Histogram of the Errors for', 'Least Square Solution with Outliers'})

%Histogram without outliers
idx=find(res_ls>12);
res_ls_wo=res_ls;
res_ls_wo(idx)=[];
n2=length(res_ls_wo);

subplot(1,2,2)
hist(res_ls_wo,50);
title({'Histogram of the Errors for','Least Square Solution without Outliers'})
set(gcf,'units','points','position',[10,10,1000,350])

%Define new variables and remove outliers from this data
ti_wo = ti;
y_obs_wo = y_obs;
y_ls_wo = y_ls;
y_true_wo = y_true;

ti_wo(idx)=[];
y_obs_wo(idx)=[];
y_ls_wo(idx)=[];
y_true_wo(idx)=[];

%Goodness of the fit

%Predictions
std_noise = sqrt(sum(res_ls_wo.^2)/(n2-p));
PI_data_ls = y_ls_wo + tinv([0.025  0.975],n2-p) * std_noise;

%Investigation additional prediction interval methods. The former yields
%the same result and the latter yields a very narrow interval. They have
%been plotted to compare.

%{
PI_data_ls_wo = y_ls_wo + tinv([0.025  0.975],n2-p-2) .* (std_noise * sqrt( 1 + 1/(n2-p) + (ti_wo-mean(ti_wo)).^2 / sum( (ti_wo-mean(ti_wo)).^2 ) ) );

for i = 1:n2
PI_data_ls_wo2(i,:) = y_ls_wo(i) + tinv([0.025 0.975],n2-p) * ( [ti_wo(i); 1]'*inv(H)*[ti_wo(i); 1] );
end

figure(5)
plot(ti_wo,y_true_wo,'r',ti_wo,y_ls_wo,'g',ti_wo,PI_data_ls(:,1),'--g',ti_wo,PI_data_ls(:,2),'--g')
figure(6)
plot(ti_wo,y_true_wo,'r',ti_wo,y_ls_wo,'g',ti_wo,PI_data_ls_wo(:,1),'--g',ti_wo,PI_data_ls_wo(:,2),'--g')
figure(7)
plot(ti_wo,y_true_wo,'r',ti_wo,y_ls_wo,'g',ti_wo,PI_data_ls_wo2(:,1),'--g',ti_wo,PI_data_ls_wo2(:,2),'--g')
%}

figure(4)
plot(ti_wo,y_obs_wo,'bo',ti_wo,y_true_wo,'r',ti_wo,y_ls_wo,'g',ti_wo,PI_data_ls(:,1),'--g',ti_wo,PI_data_ls(:,2),'--g');
title('Least Squares Solution without Outliers including Prediction Interval')
legend('Observations','True Model','Least Square Solution','Low Bound PI',' High Bound PI')
xlabel('t')
ylabel('y')


%Confidence Interval for the Parameters in x
cov_param = std_noise^2*inv(H);
std_param = diag(sqrt(cov_param));
for i = 1:2
    confint_param_ls(i,:) = x_ls(i)+tinv([0.025  0.975],p)*std_param(i);
end
disp('LS: Confidence Interval for Parameters')
disp(confint_param_ls)

T = table(x_ls, tinv(0.975,p)*std_param, std_param, corrcov(cov_param));
T.Properties.VariableNames = {'Estimate','ConfidenceInterval','Variance','CorrelationMatrix'};
T.Properties.RowNames = {'x1','x2'};
disp(T)

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
PI_data_l1 = y_l1 + tinv([0.025  0.975],n2-p) * std_noise;

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

%% Q4 L-inf Estimation



