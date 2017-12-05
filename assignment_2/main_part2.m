%% PROBLEM 2
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 2.1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all 
close all
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

A = [1./x, ones(length(x),1)];
b = 1./y;

x_star = A\b;

%Compute the coefficients by substitution into the model

theta_ls(1) = inv(x_star(2));
theta_ls(2) = x_star(1)*theta_ls(1);

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




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 2.2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA



n_obs=length(x);
cs=x;%ubstrat concentration
r=y;%reaction rate
p=2;


%% Contour Plot


theta1 = 0:0.005:0.4;
theta2 = 0:0.005:5;
[Theta1,Theta2] = meshgrid(theta1,theta2);
F=zeros(length(theta2),length(theta1));
for i=1:n_obs
    F=F+(r(i)-Theta1*cs(i)./(Theta2+cs(i))).^2;
end

F=1/2*F;

figure('DefaultAxesFontSize',16)
v = [0:0.005:0.4 0:0.005:5 0:0.0001:1];
[c,h]=contour(Theta1,Theta2,F,v,'linewidth',2);
hold on;
scatter(theta_ls(1),theta_ls(2),100,'green','filled','h')
dx = 0.01; dy =0.1;
text(theta_ls(1), theta_ls(2)+dy, '$\theta_{LS}$','Color','green','FontSize',14,'Interpreter','Latex');
title('Contour plot of $\Phi$','Interpreter','Latex')
xlabel('$\theta_1$','Interpreter','Latex')
ylabel('$\theta_2$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex')

colorbar

%% LSQNONLIN  
f_=@(theta_)(r-theta_(1)*cs/(theta_(2)+cs));

x0=[0.1,1.4];

%% %%%%%%%%%%%%%%%%%%%%%%%%% GAUSS NEWTON %%%%%%

%?LargeScale?=?off? and ?LevenbergMarquardt?=?off?
%options = optimset('LargeScale','off','levenberg-marquardt','off');
%x_opt_GN = lsqnonlin(f_,x0,[],[],options);


%% %%%%%%%%%%%%%%%%%%%%%%%%% LEVEMBERG MARQUARDT %%%%%%
%?Jacobian?=?on? and ?Algorithm?=?levenberg-marquardt? 
options=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
theta_LM = lsqnonlin(f_,x0,[],[],options);

%% Error estimation
var_est=1/(n_obs-p)*sum((r-theta_LM(1)*cs./(theta_LM(2)+cs)).^2);
J=[-cs./(theta_LM(2)+cs),+theta_LM(1)*cs./(theta_LM(2)+cs).^2];
H=J'*J;
C=inv(H);
PI_theta=theta_LM+tinv([0.025  0.975],n_obs-p)'*sqrt(var_est)*sqrt(diag(C))';
r_est=theta_LM(1)*cs./(theta_LM(2)+cs);

fprintf('Variance of Measurement noise is :%f\n',var_est)
tdist=tinv([0.025  0.975],n_obs-p);
%Display the final table containing all relevant data
T_lm = table(theta_LM',  tdist(2) *sqrt(var_est)*sqrt(diag(C)), C);
T_lm.Properties.VariableNames = {'Estimate','ConfidenceInterval','CovarianceMatrix'};
T_lm.Properties.RowNames = {'theta_1','theta_2'};
disp(T_lm)
%}

%% Final Contour Plot
figure('DefaultAxesFontSize',16)
[c,h]=contour(Theta1,Theta2,F,v,'linewidth',2);
hold on;
scatter([theta_ls(1);theta_LM(1)],[theta_ls(2);theta_LM(2)],100,'green','filled','h')
dx = -0.01; dy =0.2;
text([theta_ls(1);theta_LM(1)+dx], [theta_ls(2)+dy;theta_LM(2)+dy], ['$\theta_{LS}$';'$\theta_{LM}$'],'Color','green','FontSize',14,'Interpreter','Latex');
title('Contour plot of $\Phi$','Interpreter','Latex')
xlabel('$\theta_1$','Interpreter','Latex')
ylabel('$\theta_2$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex')
colorbar


%% Plotting the model prediction against the data

%Constructing model data from least squares estimate
x_new = 0.0:0.01:2.5;
y_model_LS = ( theta_ls(1)*x_new ) ./ (theta_ls(2)+x_new);
y_model_LM = ( theta_LM(1)*x_new ) ./ (theta_LM(2)+x_new);

figure('DefaultAxesFontSize',16)
plot(x,y,'bo',x_new,y_model_LS,'r',x_new,y_model_LM,'g','LineWidth',1.2);
title('Reaction Rate $y$ vs Substrate Concentration $x$','Interpreter','Latex')
legend({'Data','Model Prediction LS','Model Prediction LM'},'Location','northwest')
xlabel('Substrate Concentration $x$','Interpreter','Latex')
ylabel('Reaction Rate $y$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex')



%% %%%%%%%%%%%%%%%%%%%%%%%%% Best Starting Guess %%%%%%
%?Jacobian?=?on? and ?Algorithm?=?levenberg-marquardt? 
options=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
x0_b=[0.1,1.4];
[theta_LM_b,~,~,~,OUTPUT] = lsqnonlin(f_,x0_b,[],[],options);

fprintf('The number of iterations before convergence is %d with %d functions evaluated\n',...
    OUTPUT.iterations,OUTPUT.funcCount);