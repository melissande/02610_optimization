function dxdt = diffmodel(t,x,p)
%Function that computes the differential system dx(t)/dt. p is the
%parameters.

dxdt = - ( p(1) * x ) / ( p(2) + x );
end