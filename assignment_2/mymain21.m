%% Problem 2 PART 1

close all
clear all
clc

%% Initial

x = [2,2,0.667,0.667,0.4,0.4,0.286,0.286,0.222,0.220,0.2,0.2]';
y = [0.0615,0.0527,0.0334,0.0334,0.0138,0.0258,0.0129,0.0183,0.0083,0.0169,0.0129,0.0087]';

invx = 1./x;
b = 1./y;

%% Plotting (x,y) and (1/x,1/y)

figure('DefaultAxesFontSize',16)
plot(x,y,'bo');
title('Reaction Rate $y$ vs Substrate Concentration $x$','Interpreter','Latex')
legend('Data','Location','northwest')
xlabel('Substrate Concentration $x$','Interpreter','Latex')
ylabel('Reaction Rate $y$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex')
set(gcf,'units','points','position',[10,10,900,450])
xlim([0 2.05])

figure('DefaultAxesFontSize',16)
plot(invx,b,'bo')
title('Inv. Reaction Rate $\frac{1}{y}$ vs Inv. Substrate Concentration $\frac{1}{x}$','Interpreter','Latex')
xlabel('Inv. Substrate Concentration $\frac{1}{x}$','Interpreter','Latex')
ylabel('Inv. Reaction Rate $\frac{1}{y}$','Interpreter','Latex')
legend('Data','Location','northwest')
%set(gcf,'units','points','position',[10,10,500,12500])
%set(findobj(gcf,'type','axes'),'TickLabelInterpreter','Latex')
set(gcf,'units','points','position',[10,10,900,450])
set(gca,'TickLabelInterpreter','Latex')
xlim([0 5.1])

%% Linear least squares solution

%Solve the linear system of equations

A = [invx, ones(length(x),1)];
x_star = A\b;

b_ls = A*x_star;
r_ls = b - b_ls;

%Compute the coefficients by substitution into the model

theta_ls(1) = inv(x_star(2));
theta_ls(2) = x_star(1)*theta_ls(1);

%% Plotting prediction intervals for the linear model

n = length( r_ls );
p = 2;

std_noise_ls = sqrt( sum( r_ls.^2 ) / ( n - p ) );
for i = 1:n
PI(i,:) = b_ls(i) + tinv( [0.025 0.975],n-p ) * std_noise_ls *sqrt( 1 + [invx(i); 1]'*inv(A'*A)*[invx(i); 1] );
end

figure('DefaultAxesFontSize',16)
plot(invx,b,'bo',invx,b_ls,'r-',invx,PI(:,2),'g--',invx,PI(:,1),'g--','LineWidth',1.5)
title('Inv. Reaction Rate $\frac{1}{y}$ vs Inv. Substrate Concentration $\frac{1}{x}$','Interpreter','Latex')
xlabel('Inv. Substrate Concentration $\frac{1}{x}$','Interpreter','Latex')
ylabel('Inv. Reaction Rate $\frac{1}{y}$','Interpreter','Latex')
legend({'Data','Model Prediction','Lower Prediction Interval','Upper Prediction Interval'},'Location','northwest')
set(gcf,'units','points','position',[10,10,900,450])
set(findobj(gcf,'type','axes'),'TickLabelInterpreter','Latex')
xlim([0.4 5.1])

%% Plotting the model prediction against the data

%Constructing model data from least squares estimate
x_new = 0.0:0.01:2.5;
y_model = ( theta_ls(1)*x_new ) ./ (theta_ls(2)+x_new);

figure('DefaultAxesFontSize',16)
plot(x,y,'bo',x_new,y_model,'r','LineWidth',1.5);
title('Reaction Rate $y$ vs Substrate Concentration $x$','Interpreter','Latex')
legend({'Data','Model Prediction'},'Location','northwest')
xlabel('Substrate Concentration $x$','Interpreter','Latex')
ylabel('Reaction Rate $y$','Interpreter','Latex')
set(gcf,'units','points','position',[10,10,900,450])
set(gca,'TickLabelInterpreter','Latex')
