function [ xn,fn,gn,dv ] = my_line_search(fun,x ,f,g,h,opts,varargin )
%my_line_search  Find  am = argmin_{a > 0}{ P(a) = f(x+a*h) } , where  x  and  
% h  are given n-vectors and the scalar function  f  and its gradient  g  
% (with elements  g(i) = Df/Dx_i ) must be given by a MATLAB function with 
% declaration
%            function  [f, g] = fun(x,p1,p2,...)
% p1,p2,... are parameters of the function.
%
%  Call  [xn,fn,gn,dv] = linesearch(fun,x,f,g,h)
%        [xn,fn,gn,dv] = linesearch(fun,x,f,g,h,opts,p1,p2,...)
%        [xn,fn,gn,dv] = linesearch(......)
%
% Input parameters
% fun  :  Handle to the function.
% x    :  Current x.
% f,g  :  f(x) and g(x).
% h    :  Step vector.
% opts :  Either a struct with fields  'beta1', 'beta2', 
%         'amax'  or a vector with the values of these options,
%         opts = [beta1  beta2   amax].
%         beta1, beta2   :  options for stopping criteria.
%                          P(a) <= P(0) + a*beta1*P'(0)  and  
%                          P'(a) >= beta2*P'(0).  
%                          Default  beta1 = 1e-3,  beta2 = 0.99 .
%         amax       :  Maximal allowable step.   Default  amax = 10.
%         kmax       : max number of iterations. Defulat kmax=50
% p1,p2,..  are passed dirctly to the function FUN .    
%
% Output parameters
% xn    :  x + am*h
% fn,gn :  f(xn) and g(xn).
% dv: function eval

%Default assignments
% Check OPTS
if  (nargin < 5 || isempty(opts))
    amax=10;
    beta1=1e-3;
    beta2=0.99;
    kmax=50;
else
    beta1 = opts(1);  beta2 = opts(2);  amax = opts(3); kmax=opts(4); 
end



k=0;
dv=0;


if (h'*g<0)
        alpha=0;a=0;b=min(1,amax);stop=false;
        
        while (~stop && k<kmax)
            k=k+1;
            x_bh_cell=num2cell(x+b*h);
            [f_b,df_b]=feval(fun,x_bh_cell{:},varargin{:});
            if f_b<f+beta1*h'*g*b %P(b)<lambda(b)
                %lambda(alpha)<-P(0)+beta1*P'(0)*alpha
                %P(alpha)<-f(x+alpha*h)
                a=b;
                if (h'*df_b<beta2*h'*g && b<amax)%P'(b)<beta2*P'(0) and b<amax
                    b=min(2*b,amax);
                else
                    stop=true;
                end
                
            elseif (a==0 && h'*df_b<0) %a==0 and P'(b)<0
                %P'(alpha)=h^t*f(x+alpha*h)'
                b=b/10;
            else
                stop=true; 
            end    
        end
        
        alpha=b;
        stop=((a>0 &&  h'*df_b>=beta2*h'*g)||(b>=amax && h'*df_b<beta2*h'*g));
        
        while (~stop && k<kmax)
            k=k+1;
            [alpha,a,b]=refine(a,b,x,fun,f,g,h);
             x_alpha_cell=num2cell(x+alpha*h);
            [f_alpha,df_alpha]=feval(fun,x_alpha_cell{:},varargin{:});
            stop=((f_alpha<=f+beta1*h'*g*alpha)&&(h'*df_alpha>=beta2*h'*g));
            
        end
        x_alpha_cell=num2cell(x+alpha*h);
        [f_alpha,df_alpha]=feval(fun,x_alpha_cell{:},varargin{:});
        if f_alpha>=f
            alpha=0;
        end
        
        xn=x+alpha*h;
        xn_cell=num2cell(xn);
        [fn,gn]=feval(fun,xn_cell{:},varargin{:});
    
else disp('Descent direction is not downhill - END') 
end


end

