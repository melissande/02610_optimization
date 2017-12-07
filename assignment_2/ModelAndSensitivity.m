function dzdt = ModelAndSensitivity(t,z,p,n)
%Function that allows the computation of the both differential equation and
%the sensitivities of the model function. p is the parameters and n is the
%size of the x(t) system of differential equations.



%The input vector is sorted out. This is necessary to compute the
%sensitivities
x = z(1:n,1);
sp = z(n+1:end,1);

%Compuation of the derivatives necessary
dxdt = - ( p(1) * x ) / ( p(2) + x );

dfdx = - p(1) / ( p(2)+x ) + p(1)*x / ( p(2)+x )^2;
dfdp = [ -x / (p(2) +x) ; p(1)*x / ( p(2)+x )^2 ];

%The second system of differential equations for the time derivative of
%the sensitivity

Spdot = dfdx*sp + dfdp;

dzdt = [dxdt; Spdot];        %Return derivatives as a vector
