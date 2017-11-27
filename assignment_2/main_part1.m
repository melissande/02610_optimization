%% PROBLEM 1
close all
clear all
clc

%% Q1: Load Data and establish model

%Load Data

data_pb1=load('Problem1Data.mat');

%Construct model vectors etc.

ti=data_pb1.t;
y_obs=data_pb1.y;

alpha=1.0;
beta=0.0;

y_true=alpha*ti+beta;

% Initiatial Variables
n=length(y_true);
A=[ti,ones(n,1)];
b=y_obs;
p=2;
idx = y_obs > 12;

%% Plot data

figure('DefaultAxesFontSize',16)
plot(ti,y_obs,'bo',ti,y_true,'r');
title('Observations and True Model','Interpreter','Latex')
legend('Observations','True Model','Location','northwest')
xlabel('t')
ylabel('y')
set(gca,'TickLabelInterpreter','Latex')


%% Q2: l2 Squares Estimation

%Constructing and solving the quadratic unconstrained least squares
%solution to the parameters of the model
H=A'*A;
g=-A'*b;
gamma=-1/2*(b)'*b;

x_ls=-H\g;
y_ls=A*x_ls;

%Plotting the found parameter-solution against the true model
figure('DefaultAxesFontSize',16)
plot(ti,y_obs,'bo',ti,y_true,'r',ti,y_ls,'g');
title('$\ell_{2}$-Solution with Outliers','Interpreter','Latex')
legend({'Observations','True Model','$\ell_{2}$-Model Solution'},'Interpreter','Latex','Location','northwest')
xlabel('t')
ylabel('y')
set(gca,'TickLabelInterpreter','Latex')

%Histogram with outliers
res_ls=y_obs-y_ls;

figure('DefaultAxesFontSize',16)
subplot(1,2,1)
hist(res_ls,50);
title({'Histogram of the Errors for', '$\ell_{2}$-Solution with Outliers'},'Interpreter','Latex')
xlabel('Residual Value')
ylabel('Counts')

%Histogram without outliers
res_ls_wo=res_ls;
res_ls_wo(idx)=[];
n2=length(res_ls_wo);

subplot(1,2,2)
hist(res_ls_wo,50);
title({'Histogram of the Errors for','$\ell_{2}$-Solution without Outliers'},'Interpreter','Latex')
xlabel('Residual Value')
ylabel('Counts')
set(gcf,'units','points','position',[10,10,1000,350])
set(findobj(gcf,'type','axes'),'TickLabelInterpreter','Latex')

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
std_noise_ls = sqrt(sum(res_ls_wo.^2)/(n2-p));
for i = 1:n2
PI_data_ls(i,:) = y_ls_wo(i) + tinv([0.025 0.975],n2-p) * std_noise_ls *sqrt( 1 + [ti_wo(i); 1]'*inv(A'*A)*[ti_wo(i); 1] );
end

%{
PI_data_ls = y_ls_wo + tinv([0.025  0.975],n2-p) * std_noise_ls;

%Below is an investigation of additional prediction interval methods. 
%The former yields the same result and the latter yields a very narrow 
%interval. They have been plotted to compare. The initial method was chosen


PI_data_ls_wo = y_ls_wo + tinv([0.025  0.975],n2-p-2) .* (std_noise_ls * sqrt( 1 + 1/(n2-p) + (ti_wo-mean(ti_wo)).^2 / sum( (ti_wo-mean(ti_wo)).^2 ) ) );

for i = 1:n2
PI_data_ls_wo2(i,:) = y_ls_wo(i) + tinv([0.025 0.975],n2-p) * std_noise_ls *sqrt( 1 + [ti_wo(i); 1]'*inv(H)*[ti_wo(i); 1] );
end

figure(51)
plot(ti_wo,y_true_wo,'r',ti_wo,y_ls_wo,'g',ti_wo,PI_data_ls(:,1),'--g',ti_wo,PI_data_ls(:,2),'--g')
figure(52)
plot(ti_wo,y_true_wo,'r',ti_wo,y_ls_wo,'g',ti_wo,PI_data_ls_wo(:,1),'--g',ti_wo,PI_data_ls_wo(:,2),'--g')
figure(53)
plot(ti_wo,y_true_wo,'r',ti_wo,y_ls_wo,'g',ti_wo,PI_data_ls_wo2(:,1),'--g',ti_wo,PI_data_ls_wo2(:,2),'--g')
%}

figure('DefaultAxesFontSize',16)
plot(ti_wo,y_obs_wo,'bo',ti_wo,y_true_wo,'r',ti_wo,y_ls_wo,'g',ti_wo,PI_data_ls(:,1),'--g',ti_wo,PI_data_ls(:,2),'--g');
title({'$\ell_{2}$-Solution without Outliers', 'including Prediction Interval'},'Interpreter','Latex')
legend({'Observations','True Model','$\ell_{2}$-Solution','Low Bound PI',' High Bound PI'},'Interpreter','Latex','Location','northwest')
xlabel('t')
ylabel('y')
set(gca,'TickLabelInterpreter','Latex')

%Confidence Interval for the Parameters in x
cov_param_ls = std_noise_ls^2*inv(A'*A);
std_param_ls = diag(sqrt(cov_param_ls));
for i = 1:2
    confint_param_ls(i,:) = x_ls(i)+tinv([0.025  0.975],n2-p)*std_param_ls(i);
end
disp('LS: Confidence Interval for Parameters')
disp(confint_param_ls)

%Display the final table containing all relevant data
T_ls = table(x_ls, tinv(0.975,n2-p)*std_param_ls, std_param_ls, corrcov(cov_param_ls));
T_ls.Properties.VariableNames = {'Estimate','ConfidenceInterval','StandardDeviation','CorrelationMatrix'};
T_ls.Properties.RowNames = {'x1','x2'};
disp(T_ls)
%}
%% Q3: l1  Estimation

%Constructing and solving the linear unconstrained program solution
%to the parameters of the model
f_l1=[zeros(p,1);ones(n,1)];
A_l1=[-A,-eye(n);A,-eye(n)];
b_l1=[-b;b];
  
x_l1 = linprog(f_l1,A_l1,b_l1); 
x_l1=x_l1(1:2);
y_l1=A*x_l1;

%Plotting the found parameter-solution against the true model
figure('DefaultAxesFontSize',16)
plot(ti,y_obs,'bo',ti,y_true,'r',ti,y_l1,'g');
title('$\ell_{1}$-Solution with Outliers','Interpreter','Latex')
legend({'Observations','True Model','$\ell_{1}$-Model Solution'},'Interpreter','Latex','Location','northwest')
xlabel('t')
ylabel('y')
set(gca,'TickLabelInterpreter','Latex')

%Histogram with outliers
res_l1=y_obs-y_l1;
figure('DefaultAxesFontSize',16)
subplot(1,2,1)
hist(res_l1,50);
title({'Histogram of the Errors for', '$\ell_{1}$-Solution with Outliers'},'Interpreter','Latex')
xlabel('Residual Value')
ylabel('Counts')

%Histogram without outliers
res_l1_wo=res_l1;
res_l1_wo(idx)=[];
n2=length(res_l1_wo);

subplot(1,2,2)
hist(res_l1_wo,50);
title({'Histogram of the Errors for', '$\ell_{1}$-Solution without Outliers'},'Interpreter','Latex')
xlabel('Residual Value')
ylabel('Counts')
set(gcf,'units','points','position',[10,10,1000,350])
set(findobj(gcf,'type','axes'),'TickLabelInterpreter','Latex')


%Define new variables and remove outliers from this data
ti_wo = ti;
y_obs_wo = y_obs;
y_l1_wo = y_l1;
y_true_wo = y_true;

ti_wo(idx)=[];
y_obs_wo(idx)=[];
y_l1_wo(idx)=[];
y_true_wo(idx)=[];

%Goodness of the fit

%Predictions
std_noise_l1 = sqrt(sum(res_l1_wo.^2)/(n2-p));
%PI_data_l1 = y_l1_wo + tinv([0.025  0.975],n2-p) * std_noise_l1;
for i = 1:n2
PI_data_l1(i,:) = y_l1_wo(i) + tinv([0.025 0.975],n2-p) * std_noise_l1 *sqrt( 1 + [ti_wo(i); 1]'*inv(A'*A)*[ti_wo(i); 1] );
end

figure('DefaultAxesFontSize',16)
plot(ti_wo,y_obs_wo,'bo',ti_wo,y_true_wo,'r',ti_wo,y_l1_wo,'g',ti_wo,PI_data_l1(:,1),'--g',ti_wo,PI_data_l1(:,2),'--g');
title({'$\ell_{1}$-Solution without Outliers', 'including Prediction Interval'},'Interpreter','Latex')
legend({'Observations','True Model','$\ell_{1}$-Solution','Low Bound PI',' High Bound PI'},'Interpreter','Latex','Location','northwest')
xlabel('t')
ylabel('y')
set(gca,'TickLabelInterpreter','Latex')

%Confidence Interval for the Parameters in x
cov_param_l1 = std_noise_l1^2*inv(A'*A);
std_param_l1 = diag(sqrt(cov_param_l1));
for i = 1:2
    confint_param_l1(i,:) = x_l1(i)+tinv([0.025  0.975],n2-p)*std_param_l1(i);
end
disp('LS: Confidence Interval for Parameters')
disp(confint_param_l1)

%Display the final table containing all relevant data
T_l1 = table(x_l1, tinv(0.975,n2-p)*std_param_l1, std_param_l1, corrcov(cov_param_l1));
T_l1.Properties.VariableNames = {'Estimate','ConfidenceInterval','StandardDeviation','CorrelationMatrix'};
T_l1.Properties.RowNames = {'x1','x2'};
disp(T_l1)

%% Q4 L-inf Estimation

%Constructing and solving the linear unconstrained program solution
%to the parameters of the model
f_linf = [zeros(p,1);1];
A_linf = [-A,-ones(n,1);A,-ones(n,1)];
b_linf = [-b;b];

x_linf = linprog(f_linf,A_linf,b_linf);
x_linf = x_linf(1:2);
y_linf = A*x_linf;

%Plotting the found parameter-solution against the true model
figure('DefaultAxesFontSize',16)
plot(ti,y_obs,'bo',ti,y_true,'r',ti,y_linf,'g');
title('$L_{\infty}$-Solution with Outliers','Interpreter','Latex')
legend({'Observations','True Model','$\ell_{\infty}$-Model Solution'},'Interpreter','Latex','Location','northwest')
xlabel('t')
ylabel('y')
set(gca,'TickLabelInterpreter','Latex')

%Histogram with outliers
res_linf = y_obs-y_linf;

figure('DefaultAxesFontSize',16)
subplot(1,2,1)
hist(res_linf,50);
title({'Histogram of the Errors for', '$\ell_{\infty}$-Solution with Outliers'},'Interpreter','Latex')
xlabel('Residual Value')
ylabel('Counts')

%Histogram without outliers
res_linf_wo = res_linf;
res_linf_wo(idx)=[];
n2 =length(res_linf_wo);

subplot(1,2,2)
hist(res_linf_wo,50);
title({'Histogram of the Errors for','$\ell_{\infty}$-Solution without Outliers'},'Interpreter','Latex')
xlabel('Residual Value')
ylabel('Counts')
set(gcf,'units','points','position',[10,10,1000,350])
set(findobj(gcf,'type','axes'),'TickLabelInterpreter','Latex')

%Define new variables and remove outliers from this data
ti_wo = ti;
y_obs_wo = y_obs;
y_linf_wo = y_linf;
y_true_wo = y_true;

ti_wo(idx)=[];
y_obs_wo(idx)=[];
y_linf_wo(idx)=[];
y_true_wo(idx)=[];

%Goodness of the fit

%Predictions
std_noise_linf = sqrt(sum(res_linf_wo.^2)/(n2-p));
for i = 1:n2
PI_data_linf(i,:) = y_linf_wo(i) + tinv([0.025 0.975],n2-p) * std_noise_linf *sqrt( 1 + [ti_wo(i); 1]'*inv(A'*A)*[ti_wo(i); 1] );
end
%PI_data_linf = y_linf_wo + tinv([0.025  0.975],n2-p) * std_noise_linf;

figure('DefaultAxesFontSize',16)
plot(ti_wo,y_obs_wo,'bo',ti_wo,y_true_wo,'r',ti_wo,y_linf_wo,'g',ti_wo,PI_data_linf(:,1),'--g',ti_wo,PI_data_linf(:,2),'--g');
title({'$\ell_{\infty}$-Solution without Outliers', 'including Prediction Interval'},'Interpreter','Latex')
legend({'Observations','True Model','$\ell_{\infty}$-Solution','Low Bound PI',' High Bound PI'},'Interpreter','Latex','Location','northwest')
xlabel('t')
ylabel('y')
set(gca,'TickLabelInterpreter','Latex')

%Confidence Interval for the Parameters in x
cov_param_linf = std_noise_linf^2*inv(A'*A);
std_param_linf = diag(sqrt(cov_param_linf));
for i = 1:2
    confint_param_linf(i,:) = x_linf(i)+tinv([0.025  0.975],n2-p)*std_param_linf(i);
end
disp('LS: Confidence Interval for Parameters')
disp(confint_param_linf)

%Display the final table containing all relevant data
T_linf = table(x_linf, tinv(0.975,n2-p)*std_param_linf, std_param_linf, corrcov(cov_param_linf));
T_linf.Properties.VariableNames = {'Estimate','ConfidenceInterval','StandardDeviation','CorrelationMatrix'};
T_linf.Properties.RowNames = {'x1','x2'};
disp(T_linf)

%% Huber-Estimation

%Constructing and solving the constrained quadratic program solution
%to the parameters of the model

%Refine A to contain 100 columns of zeros.
A_new = [A,zeros(100,100-p)];
tau = 3;

%Define necessary inputs to quadprog
H_huber = [zeros(3*n,n),zeros(3*n,n),zeros(3*n,n),zeros(3*n,n);zeros(n,n),zeros(n,n),zeros(n,n),eye(n)];
f_huber = [zeros(n,1);tau*ones(n,1);tau*ones(n,1);zeros(n,1)];
A_huber = [-A_new,-eye(n),eye(n),eye(n)];
b_huber = -b;
LB = [-Inf(n,1);zeros(n,1);zeros(n,1);-Inf(n,1)];
UB = Inf(4*n,1);

x_huber = quadprog(H_huber,f_huber,A_huber,b_huber,A_huber,b_huber,LB,UB);
x_huber = x_huber(1:2);
y_huber= A*x_huber;

%Plotting the found parameter-solution against the true model
figure('DefaultAxesFontSize',16)
plot(ti,y_obs,'bo',ti,y_true,'r',ti,y_huber,'g');
title('Huber-Solution with Outliers','Interpreter','Latex')
legend({'Observations','True Model','Huber-Model Solution'},'Interpreter','Latex','Location','northwest')
xlabel('t')
ylabel('y')
set(gca,'TickLabelInterpreter','Latex')


%Histogram with outliers
res_huber = y_obs-y_huber;

figure('DefaultAxesFontSize',16)
subplot(1,2,1)
hist(res_huber,50);
title({'Histogram of the Errors for', 'Huber-Solution with Outliers'},'Interpreter','Latex')
xlabel('Residual Value')
ylabel('Counts')

%Histogram without outliers
res_huber_wo = res_huber;
res_huber_wo(idx)=[];
n2 =length(res_huber_wo);

subplot(1,2,2)
hist(res_huber_wo,50);
title({'Histogram of the Errors for','Huber-Solution without Outliers'},'Interpreter','Latex')
xlabel('Residual Value')
ylabel('Counts')
set(gcf,'units','points','position',[10,10,1000,350])
set(findobj(gcf,'type','axes'),'TickLabelInterpreter','Latex')

%Define new variables and remove outliers from this data
ti_wo = ti;
y_obs_wo = y_obs;
y_huber_wo = y_huber;
y_true_wo = y_true;

ti_wo(idx)=[];
y_obs_wo(idx)=[];
y_huber_wo(idx)=[];
y_true_wo(idx)=[];

%Goodness of the fit

%Predictions
std_noise_huber = sqrt(sum(res_huber_wo.^2)/(n2-p));
for i = 1:n2
PI_data_huber(i,:) = y_huber_wo(i) + tinv([0.025 0.975],n2-p) * std_noise_huber *sqrt( 1 + [ti_wo(i); 1]'*inv(A'*A)*[ti_wo(i); 1] );
end
%PI_data_huber = y_huber_wo + tinv([0.025  0.975],n2-p) * std_noise_huber;

figure('DefaultAxesFontSize',16)
plot(ti_wo,y_obs_wo,'bo',ti_wo,y_true_wo,'r',ti_wo,y_huber_wo,'g',ti_wo,PI_data_huber(:,1),'--g',ti_wo,PI_data_huber(:,2),'--g');
title({'Huber-Solution without Outliers', 'including Prediction Interval'},'Interpreter','Latex')
legend({'Observations','True Model','Huber-Solution','Low Bound PI',' High Bound PI'},'Interpreter','Latex','Location','northwest')
xlabel('t')
ylabel('y')
set(gca,'TickLabelInterpreter','Latex')

%Confidence Interval for the Parameters in x
cov_param_huber = std_noise_huber^2*inv(A'*A);
std_param_huber = diag(sqrt(cov_param_huber));
for i = 1:2
    confint_param_huber(i,:) = x_huber(i)+tinv([0.025  0.975],n2-p)*std_param_huber(i);
end
disp('LS: Confidence Interval for Parameters')
disp(confint_param_huber)

%Display the final table containing all relevant data
T_huber = table(x_huber, tinv(0.975,n2-p)*std_param_huber, std_param_huber, corrcov(cov_param_huber));
T_huber.Properties.VariableNames = {'Estimate','ConfidenceInterval','StandardDeviation','CorrelationMatrix'};
T_huber.Properties.RowNames = {'x1','x2'};
disp(T_huber)
