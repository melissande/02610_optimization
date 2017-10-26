function [x_star,r_star] = NOfit1(t,y,n)
%% Description
%NOfit function solves the linear least squares problem in order to
%estimate the best model fitting-parameters to the data points presented in
%the vectors t and y. The order of the fit is given by n, which must be an
%odd number. The underlying basis functions are the trigonometric functions

%% INPUT
%t: time in hours
%y: measurements of NO concentration in ug/m3 
%n: order of the model

%% OUTPUT
%x_star: vector of best-fit parameters
%r_star: residual vector given the parameters in x_star

%%
%Initial Parameters
w = 2*pi/24;

%Defining the A matrix
A = zeros(size(y,1),n);

%First Column of 1's
A(:,1) = 1;

%Subsequently Columns: Trigonometric Functions
for i=1:floor(n/2)
    A(:,2*i:2*i+1)=[sin(i*w*t),cos(i*w*t)];
end

%Solving the Linear System of Equations

x_star = A\y;
r_star = y-A*x_star;

end

