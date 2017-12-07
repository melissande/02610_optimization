function [r,J] = odefun(p,y)

[~,Z] = ode45(@ModelAndSensitivity,0:10:200,[10,0,0],[],p,1);

%Residual vector
y_hat = Z(:,1);
r = y_hat - y;

%Jacobian
J = Z(:,2:3);