function [f,df,J,r] = func_basis_3(x)
%% Description
%func_basis_3 returns the Function value f(x), its jacobian J(x), 
%and its gradient df(x), for the function defined in Part 2 of Assignment 1

%% INPUT
% x: Vector where f is to be evaluated (vector of two points)
%% OUTPUT 
% f: Function value
% r: Residual vector
% J: Jacobian value
% df:Gradient value

x1 = x(1,1);
x2 = x(2,1);

f = (x1^2+x2-11)^2+(x1+x2^2-7)^2;

df = [4*x1*(x1^2+x2-11)+2*(x1+x2^2-7) ; 2*(x1^2+x2-11)+4*x2*(x1+x2^2-7)];

J = sqrt(2) * [2*x1,1 ; 1,2*x2];

r = sqrt(2) * [x1^2+x2-11;x1+x2^2-7];
