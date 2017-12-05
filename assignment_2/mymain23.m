%% Problem 2.3
clear all;close all;clc

load('MMBatchData.mat')
t = data(:,1);
x = data(:,2);

%% Plotting the data

figure('DefaultAxesFontSize',16)
plot(t,x,'bo');
title('Substrate Concentration $x(t)$ as function of time $t$','Interpreter','Latex')
legend('Data','Location','northeast')
xlabel('Time $t$','Interpreter','Latex')
ylabel('Substrate Concentration $x(t)$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex')
xlim([-5 205])
ylim([-0.5 10.5])

%% 
%Solve the differential equation using ode45 and calling the diff-equation. 
%The last input is the parameter values p of @diffmodel

[T,X] = ode45(@diffmodel,[0:10:200],10,[],[3,5]);

%% Contour Plot

%We feed various inputs of p to the above solution to the ode to generate
%the various y(theta,t). Phi can then be computed



