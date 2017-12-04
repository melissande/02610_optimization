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

[T,X] = ode45(@diffmodel,[0:10:200],10,[],[3,5])

