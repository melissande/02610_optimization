function dxdt = diffmodel(t,x,p)
dxdt = - ( p(2) * x ) / (p(1) + x );
end