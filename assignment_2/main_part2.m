%% PROBLEM 2
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 2.1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Q1: PLOT DATA

data_pb2=[2.0000 ,   0.0615;
    2.0000 ,   0.0527;
    0.6670 ,   0.0334;
    0.6670 ,   0.0334;
    0.4000 ,   0.0138;
    0.4000 ,   0.0258;
    0.2860 ,   0.0129;
    0.2860 ,   0.0183;
    0.2220 ,   0.0083;
    0.2200 ,   0.0169;
    0.2000 ,   0.0129;
    0.2000 ,   0.0087];

n_obs=length(data_pb2);
cs=data_pb2(:,1);%ubstrat concentration
r=data_pb2(:,2);%reaction rate
p=2;
%% Plot data and transformed data

figure(1);
subplot(1,2,1)
plot(cs,r,'bo')
xlabel('Substrate Concentration Cs')
ylabel('Reaction Rate r')
title('r=f(Cs,t)')
subplot(1,2,2)
plot(1./cs,1./r,'bo')
xlabel('1/Cs')
ylabel('1/r')
title('1/r=f(1/Cs,t)')

%% l1 estimation

%p1=theta2/theta1 and p2=1/theta1 ---> 1/r=p1*1/cs+p2

r_inv=1./r;
cs_inv=1./cs;
n=length(r_inv);
A=[cs_inv,ones(n,1)];
b=r_inv;


f_lp=[zeros(p,1);ones(n,1)];
A_lp=[-A,-eye(n);A,-eye(n)];
b_lp=[-b;b];

x_l1 = linprog(f_lp,A_lp,b_lp); 
x_l1=x_l1(1:2);
y_l1=A*x_l1;

%theta2=p1/p2 and theta1=1/p2
theta1=1/x_l1(2)
theta2=x_l1(1)/x_l1(2)

figure(2)
plot(cs_inv,r_inv,'bo',cs_inv,y_l1,'g');
legend('Observations','L1 Estimate')
title('L1 Estimate')
xlabel('1/cs')
ylabel('1/r')

res_l1=r_inv-y_l1;

%Goodness of the fit

%Predictions
std_noise=sqrt(sum(res_l1.^2)/length(res_l1));
PI_data_l1=y_l1+tinv([0.025  0.975],n-p)*std_noise;

figure(3)
plot(cs_inv,r_inv,'bo',cs_inv,y_l1,'g',cs_inv,PI_data_l1(:,1),'--g',cs_inv,PI_data_l1(:,2),'--g')
legend('Observations','L1 Estimate','Low Bound PI',' High Bound PI')
xlabel('1/cs')
ylabel('1/r')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 2.2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Contour Plot

theta1 = 0:0.005:0.4;
theta2 = 0:0.005:5;
[Theta1,Theta2] = meshgrid(theta1,theta2);
F=zeros(length(theta2),length(theta1));
for i=1:n_obs
    F=F+(r(i)-Theta1*cs(i)./(Theta2+cs(i))).^2;
end

F=1/2*sqrt(F);

 figure(4);
v = [0:0.01:5 0:0.01:0.8 0:0.005:2];
[c,h]=contour(Theta1,Theta2,F,v,'linewidth',2);
title('Contour plot of Phi')
xlabel('Theta 1')
ylabel('Theta 2')
colorbar


%% LSQNONLIN  
f_=@(theta_)(r-theta_(1)*cs/(theta_(2)+cs));

x0=[0.1,1.4];
x_opt = lsqnonlin(f_,x0);
%% Error estimation
var_est=1/(n-p)*sum((r-x_opt(1)*cs./(x_opt(2)+cs)).^2);
J=[cs./(x_opt(2)+cs),-x_opt(1)*cs./(x_opt(2)+cs).^2];
H=J'*J;

C=inv(H);
PI_theta=x_opt+[1,-1]'.* tinv([0.025  0.975],n-p) *sqrt(var_est)*sqrt(C);
r_est=x_opt(1)*cs./(x_opt(2)+cs);

%% Final Contour Plot
theta1 = 0:0.005:0.4;
theta2 = 0:0.005:5;
[Theta1,Theta2] = meshgrid(theta1,theta2);
F=zeros(length(theta2),length(theta1));
for i=1:n
 
    F=F+(r(i)-Theta1*cs(i)./(Theta2+cs(i))).^2;
end

F=1/2*sqrt(F);

 figure(5);
v = [0:0.01:5 0:0.01:0.8 0:0.005:2];
[c,h]=contour(Theta1,Theta2,F,v,'linewidth',2);
hold on;
 scatter(x_opt(1),x_opt(2),100,'reds','filled','h')
title('Contour plot of Phi')
xlabel('Theta 1')
ylabel('Theta 2')
colorbar

%% Plot reaction

