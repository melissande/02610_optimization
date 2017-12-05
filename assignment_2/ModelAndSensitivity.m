function zdot = ModelAndSensitivity(t,z,p,n,np)

x = z(1:n,1);
sp = z(n+1:end,1);
Sp = reshape(sp,n,np);

xdot = - ( p(1) * x ) / (p(2) + x );

"[dfdx,dfdp] = Derivatives(t,x,p);"  %Evaluate the derivatives

Spdot = dfdx*Sp + dfdp;

zdot = [xdot; Spdot(:)];             %Return derivatives as a vector