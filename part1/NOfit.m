function [x_star,r_star,A] = NOfit(t,y,n)
%NOfit that computes the least squares fit with the model defines in the function
%to the data points in the vectors t and y. The order of the fit is n, 
%which must be an odd number, and the underlying basis functions are the trigonometric functions
%% INPUT:
%t:time in hours
%y: measurements of NO concentration in ?g/m3 
%n:order of the model
%% OUTPUT
%x_star: model vector
%r_star:residual 
%A: matrix model
%%

A=zeros(size(y,1),n);
A(:,1)=ones(size(y,1),1);
w=2*pi/size(y,1);
for i=1:floor(n/2)
    A(:,2*i:2*i+1)=[sin(i*w*t),cos(i*w*t)];
end

x_star=A\y;
r_star=y-A*x_star;



end

