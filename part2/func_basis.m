function [ f,df ] = func_basis(x)
%func_basis returns the function,its gradient 
%for the function defined in PART 2 assignment 1
%% INPUT
%x: vector where f is to be evaluated (vector of two points)
%% OUTPUT 
%f:function
%df:gradient
%%
x1=x(1,1);
x2=x(2,1);
f =(x1^2+x2-11)^2+(x1+x2^2-7)^2;
df=[4*x1^3+4*x1*x2-42*x1+2*x2^2-14;4*x2^3+2*x1^2-26*x2+4*x1*x2-22];



end
