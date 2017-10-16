function [alpha,a,b]=refine(a,b,x,fun,f,g,h)
 %The input is an interval [a, b] which we know contains acceptable points, and
%the output is an alpha found by interpolation. We want to be sure that the
%intervals have strictly decreasing widths, so we only accept the new alpha
%if it is inside [a+d, b-d], where d = 1/10 (b-a). The alpha splits [a, b] into 10
%two subintervals, and we also return the subinterval which must contain acceptable points.
x_ah_cell=num2cell(x+a*h);
            [f_a,df_a]=feval(fun,x_ah_cell{:},varargin{:});

x_bh_cell=num2cell(x+b*h);
            [f_b,df_b]=feval(fun,x_bh_cell{:},varargin{:});

D=b-a;
c=(f_b-f_a-D*h'*df_a)/D^2; %obtained because of taylor dvp second order of phi(t)=P(a)+P'(a)*(t-a)+c*(t-a)^2
%if c>0, phi has minimizer otherwise alpha midpoint of [a,b]
if c>0
    alpha=(a-h'*df_a)/(2*c);
    alpha=min(max(alpha,a +0.1*D),b-0.1*D);
else
    alpha=(a+b)/2;
end
x_alpha_cell=num2cell(x+alpha*h);
            [f_alpha,df_alpha]=feval(fun,x_alpha_cell{:},varargin{:});
            
if f_alpha<f+beta1*h'*g*alpha
    a=alpha;
else
    b=alpha;

end

