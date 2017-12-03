%% Problem 2 PART 1

close all
clear all
clc

%% Initial

x = [2,2,0.667,0.667,0.4,0.4,0.286,0.286,0.222,0.220,0.2,0.2]';
y = [0.0615,0.0527,0.0334,0.0334,0.0138,0.0258,0.0129,0.0183,0.0083,0.0169,0.0129,0.0087]';

%% Plotting (x,y) and (1/x,1/y)

figure('DefaultAxesFontSize',16)
subplot(1,2,1)
plot(x,y,'bo');
title('Reaction Rate $y$ vs Substrate Concentration $x$','Interpreter','Latex')
legend('Data','Location','northwest')
xlabel('Substrate Concentration $x$','Interpreter','Latex')
ylabel('Reaction Rate $y$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex')

subplot(1,2,2)
plot(1./x,1./y,'bo')
title('Inv. Reaction Rate $\frac{1}{y}$ vs Inv. Substrate Concentration $\frac{1}{x}$','Interpreter','Latex')
xlabel('Inv. Substrate Concentration $\frac{1}{x}$','Interpreter','Latex')
ylabel('Inv. Reaction Rate $\frac{1}{y}$','Interpreter','Latex')
legend('Data','Location','northwest')
set(gcf,'units','points','position',[10,10,1250,350])
set(findobj(gcf,'type','axes'),'TickLabelInterpreter','Latex')

%% Linear least squares solution

%Solve the linear system of equations

A = [1./x, ones(length(x),1)];
b = 1./y;

x_star = A\b;

%Compute the coefficients by substitution into the model

theta_ls(1) = inv(x_star(1));
theta_ls(2) = x_star(2)*theta_ls(1);

%% Plotting the model prediction against the data

%Constructing model data from least squares estimate
x_new = 0.0:0.01:2.5;
y_model = ( theta_ls(1)*x_new ) ./ (theta_ls(2)+x_new);

figure('DefaultAxesFontSize',16)
plot(x,y,'bo',x_new,y_model,'r','LineWidth',1.2);
title('Reaction Rate $y$ vs Substrate Concentration $x$','Interpreter','Latex')
legend({'Data','Model Prediction'},'Location','northwest')
xlabel('Substrate Concentration $x$','Interpreter','Latex')
ylabel('Reaction Rate $y$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex')






