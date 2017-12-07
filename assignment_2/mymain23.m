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


%% Calculating phi(theta)

%Define coordinate grid to act as parameter input space.
theta1 = linspace( 0.01,1,100 );
theta2 = linspace( 0.01,3,100 );

[T1,T2] = meshgrid (theta1 , theta2 );

[rT,cT] = size( T1 );


%We feed various inputs of p to the ode-solution to generate
%the various y_hat(theta,t) outputs.

tic
for i = 1:rT
    
    for j = 1:cT
    p = [ T1(i,j) , T2(i,j) ];
    
    %Method 1. Using diffmodel
    [~,A{i,j}] = ode45( @diffmodel , 0:10:200 , 10 , [] , p );
    Y(i,j) = 0.5 *  sum( ( y-A{i,j} ).^2 );
    
    %Method 2. Using ModelAndSensitivity
    [~,K] = ode45( @ModelAndSensitivity , 0:10:200 , [10;0;0] , [] , p , 1);
    B{i,j} = K(:,1);
    Z(i,j) = 0.5 *  sum( ( y-B{i,j} ).^2 );
    
    end
    
end
toc

%%
%Create the contour plot

figure('DefaultAxesFontSize',16)
contourf(T1,T2,Y,100)
set(gcf,'units','points','position',[10,10,900,450])

figure('DefaultAxesFontSize',16)
contourf(T1,T2,Z,100)
set(gcf,'units','points','position',[10,10,900,450])


%% Solving the unconstrained optimization problem using lsqnonlin

