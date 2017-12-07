function [r,J] = odefun(p,y)

%Solve the differential equations
[~,Z] = ode45(@ModelAndSensitivity,0:10:200,[10,0,0],[],p,1);

%Extract solution and compute the residual vector
y_hat = Z(:,1);
r = y_hat - y;

%Extract Jacobian
J = Z(:,2:3);