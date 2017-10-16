function [ f,df,d2f ] = func_basis( x ,mu)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

f = x-mu*log(x);
df=1-mu./x;
d2f=mu./x.^2;


end

