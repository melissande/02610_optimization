function [ f,df ] = func_basis_2( x ,mu)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

f = x-mu*log(x);
df=1-mu./x;



end