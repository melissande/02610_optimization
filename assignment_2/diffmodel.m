function dxdt = diffmodel(t,x,p)
dxdt = - ( p(1) * x ) / ( p(2) + x );
end