function [p,c]= cvg_rate( error,fig_number )
%cvg_rate computes the convergence rate of an algorithm provided the prediction
%error vector ordered as increasing iteration
%Depending on the type of convergence (linear, superlinear, quadratic), 
%a figure is plotted to demonstrate the convergence
%fig_number: count for figure so that it doesn't overwrite
%% Define the prediction error "shifted"
ek=error;
ek(end)=[];
ek_next=error;
ek_next(1)=[];

%% Solve list square problem on cvg rate problem

b=log(ek_next);
A=[ones(length(ek),1),log(ek)];
x=A\b;
p=round(x(2),1);
c=round(exp(x(1)),1);
%% Type of convergence
iter=1:length(ek);

figure(fig_number);
if (p==1.0)
    %Linear cvg
    x=iter';
    y=ek_next./ek;
    plot(x,y)
    xlabel('iteration k')
    ylabel('e_{k+1}/e_k')
    title('Linear Convergence')
elseif (p==2.0)
    x=iter';
    y=ek_next./ek.^2;
    plot(x,y)
    xlabel('iteration k')
    ylabel('e_{k+1}/e_k^2')
    title('Quadratic Convergence')
elseif (1<p && p<2)
    x=iter';
    y=ek_next./ek;
    plot(x,y)
     xlabel('iteration k')
    ylabel('e_{k+1}/e_k')
    title('Superlinear Convergence')

else disp('Convergence not defined') 


end

end

