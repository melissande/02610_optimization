%% Problem 2.3
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


%% 
%Solve the differential equation using ode45 and calling the diff-equation. 
%The last input is the parameter values p of @diffmodel

[~,X] = ode45(@diffmodel,0:10:200,10,[],[0.5,0.5])

%% Contour Plot

%Define coordinate grid to act as parameter input space.
theta1 = linspace(0.01,1,50);
theta2 = linspace(0.01,3,50);
[T1,T2] = meshgrid(theta1,theta2);

[rT,cT] = size(T1);

%We feed various inputs of p to the ode-solution to generate
%the various y_hat(theta,t) outputs.
tic
for i = 1:rT
    for j = 1:cT
    [~,A{i,j}] = ode45(@diffmodel,0:10:200,10,[],[T1(i,j),T2(i,j)]);
    Y(i,j) =0.5*norm(y-A{i,j})^2;
    end
end
toc

%%
%Create the contour plot
%v = [0:100:0.2 100:500 0.5 500:1000];

figure('DefaultAxesFontSize',16)
contour(T1,T2,Y,100)
set(gcf,'units','points','position',[10,10,900,450])
