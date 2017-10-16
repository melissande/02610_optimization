function p= cvg_rate( error )
%cvg_rate computes the convergence rate of an algorithm provided the prediction
%error vector ordered as increasing iteration
%Depending on the type of convergence (linear, superlinear, quadratic), 
%a figure is plotted to demonstrate the convergence

%% Define the prediction error "shifted"
ek=error;
ek(end)=[];
ek_next=error;
ek_next(1)=[];

%% Solve list square problem on cvg rate problem

b=log(ek_next);
A=[ones(length(ek),1),log(ek)];
x=A\b;
p=x(2);

%% Type of convergence
iter=1:length(ek);
bool_disp=true;
figure;
if (p==1)
    %Linear cvg
    x=iter';
    y=ek_next./ek;
    plot(x,y,'.','MarkerSize',20)
    xlabel('iteration k')
    ylabel('e_{k+1}/e_k')
    title('Linear Convergence')
elseif (p==2)
    x=iter';
    y=ek_next./ek.^2;
    plot(x,y,'.','MarkerSize',20)
    xlabel('iteration k')
    ylabel('e_{k+1}/e_k^2')
    title('Quadratic Convergence')
elseif (1<p && p<2)
    x=iter';
    y=ek_next./ek;
    plot(x,y,'.','MarkerSize',20)
     xlabel('iteration k')
    ylabel('e_{k+1}/e_k')
    title('Superlinear Convergence')

else disp('Convergence not defined') 


end

end

