function [ f,df,d2f ] = func_basis_2(x)
%func_basis returns the function,its gradient and its Hessian
%for the function defined in PART 2 assignment 1
%% INPUT
%x: vector where f is to be evaluated (vector of two points)
%% OUTPUT 
%f:function
%df:gradient
%d2f:Hessian
%%
x1=x(1,1);
x2=x(2,1);
f =(x1^2+x2-11)^2+(x1+x2^2-7)^2;
df=[4*x1^3+4*x1*x2-42*x1+2*x2^2-14;4*x2^3+2*x1^2-26*x2+4*x1*x2-22];
d2f=[4*(x1.^2+x2-11)+8*x1.^2+2 4*x1+4*x2;4*x1+4*x2 4*(x1+x2.^2-7)+8*x2.^2+2];

end