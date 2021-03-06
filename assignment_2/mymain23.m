%% Problem 2.3
%% Initialization
clear all;close all;clc

load('MMBatchData.mat')
t = data(:,1);
y = data(:,2);

%% Plotting the data

figure('DefaultAxesFontSize',16)
plot(t,y,'bo');
title('Substrate Concentration $x(t)$ as function of time $t$','Interpreter','Latex')
legend('Data','Location','northeast')
xlabel('Time $t$','Interpreter','Latex')
ylabel('Substrate Concentration $x(t)$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex')
set(gcf,'units','points','position',[10,10,900,450])
xlim([-5 205])
ylim([-0.5 10.5])


%% Calculating phi(theta)

%Define coordinate grid to act as parameter input space.
theta1 = linspace( -0.1108,0.2892,100 );
theta2 = linspace( 0.2699,1.2699,100 );

[T1,T2] = meshgrid (theta1 , theta2 );
[rT,cT] = size( T1 );


%The various inputs of p are fed to the odesolver to generate
%the corresponding y_hat(theta,t) outputs.

tic
for i = 1:rT
    
    for j = 1:cT
        
    p = [ T1(i,j) , T2(i,j) ];
    
    %Method 1. Using diffmodel
    [~,A{i,j}] = ode45( @diffmodel , 0:10:200 , 10 , [] , p );
    Y(i,j) = 0.5 *  sum( ( y-A{i,j} ).^2 );
    
    end
    
end
toc

%% Creating contour plot

figure('DefaultAxesFontSize',16)
contourf( T1,T2,Y,100 )
hold on
scatter(0.0892,0.7699,'r*')
text(0.0892+0.005,0.7699+0.01, '$\theta_{lsqnonlin}$','Color','red','FontSize',18,'Interpreter','Latex');
set(gcf,'units','points','position',[10,10,900,450])
title('Contours of $\phi(\theta)$','Interpreter','Latex')
xlabel('$\theta_{1}$','Interpreter','Latex')
ylabel('$\theta_{2}$','Interpreter','Latex')
colorbar


%% Solving the unconstrained optimization problem using lsqnonlin

optim_LM = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','SpecifyObjectiveGradient',true,'FunctionTolerance',1e-9,'StepTolerance',1e-9);
optim_TR = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','SpecifyObjectiveGradient',true,'FunctionTolerance',1e-9,'StepTolerance',1e-9);

p0 = [0.1,0.1];
[theta_LM,resnorm_LM,residual_LM,~,overview_LM] = lsqnonlin( @odefun,p0,[],[],optim_LM,y);
[theta_TR,resnorm_TR,residual_TR,~,overview_TR] = lsqnonlin( @odefun,p0,[],[],optim_TR,y);

%% Statistic Calculations

%Retrieve residual and Jacobian for the minimizer
[r,J] = feval(@odefun,theta_TR,y);

%Calculate the variance estimator
m = 21;
np = 2;

noise_estimate = (1 / ( m-np )) * sum( r.^2 );

%Calculating the covariance matrix

cov_param = noise_estimate * inv( J'*J );

%Calculating the parameter confidence intervals

theta1_confint = tinv(0.975,m-np) * sqrt( noise_estimate ) * cov_param(1,1);
theta2_confint = tinv(0.975,m-np) * sqrt( noise_estimate ) * cov_param(2,2);

%Construct table with a summary of the data

T_TR = table(theta_TR', [theta1_confint;theta2_confint], [cov_param(1,1);cov_param(2,2)], cov_param);
T_TR.Properties.VariableNames = {'Estimate','ConfidenceInterval','StandardDeviation','CovarianceMatrix'};
T_TR.Properties.RowNames = {'theta1','theta2'};
disp(T_TR)
