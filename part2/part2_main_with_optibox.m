%% PART 2  

%% Function definition

func=@(x1,x2)((x1.^2+x2-11).^2+(x1+x2.^2-7).^2);

%% Q2.1: Contour Plots

%Definition of stationary points (found analytically)
x_sad=[ 0.08667,2.88425;-3.07302,-0.08135; 3.38515,0.07385;-0.12796, -1.95371];
x_min=[ -2.80511,3.13131;3.58442,-1.84812;-3.77931,-3.28318];
x_max=[-0.27084,-0.92303];
x_global_min=[3,2];

make_contour_plot_init( func, [ 0.08667,2.88425;-3.07302,-0.08135; 3.38515,0.07385;-0.12796, -1.95371],[ -2.80511,3.13131;3.58442,-1.84812;-3.77931,-3.28318],[-0.27084,-0.92303], [3,2] )

%% Q2.2: Gradient and Hessian

%% Verify that df=0 for all sationary points
for i=1:size(x_sad,1)
    [f,df,d2f]=func_basis_2(x_sad(i,:)');
   fprintf('Saddle point %d: df=[%d,%d]  and d2f=[%d,%d;%d,%d] \n',i,df,d2f)
   [V,D] = eig(d2f);
   fprintf('The eigen values are: [%d,%d]so the matrix is Indefinite\n',diag(D));
end
fprintf('\n');

for i=1:size(x_min,1)
    [f,df,d2f]=func_basis_2(x_min(i,:)');
   fprintf('Local minimum %d: df=[%d,%d]  and d2f=[%d,%d;%d,%d] \n',i,df,d2f)
   [V,D] = eig(d2f);
   fprintf('The eigen values are: [%d,%d]so the matrix is Positive Definite\n',diag(D));
end
fprintf('\n');
for i=1:size(x_max,1)
    [f,df,d2f]=func_basis_2(x_max(i,:)');
   fprintf('Local maximum %d: df=[%d,%d]  and d2f=[%d,%d;%d,%d] \n',i,df,d2f)
   [V,D] = eig(d2f);
    fprintf('The eigen values are: [%d,%d]so the matrix is Negative Definite\n',diag(D));
end
fprintf('\n');
[f,df,d2f]=func_basis_2(x_global_min');
fprintf('Global minimum %d: df=[%d,%d]  and d2f=[%d,%d;%d,%d] \n',i,df,d2f)
[V,D] = eig(d2f);
 fprintf('The eigen values are: [%d,%d]so the matrix is Negative Definite\n',diag(D));
fprintf('\n');
%% Q2.4:  Steepest Descent algorithm
fignumber=1;
%Remark: basically falls into all the stationary points ever (not max though)
%when close to but alwayse cvg to something
x_global_min=[3,2];

disp('Pretty Close to the global minimum and far away from any other stationary point')
x0=[5,5]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=steepest_optim_ls(@func_basis,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', xopt' ,fignumber )
fignumber=fignumber+1;
%%
disp('Close to Global Minimum and surrounded by stationary points')
x0=[0,0]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=steepest_optim_ls(@func_basis,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', xopt' ,fignumber )
fignumber=fignumber+1;
%%
disp('Close from Global Minimum and with some stationary points around')
x0=[-5,0]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=steepest_optim_ls(@func_basis,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', xopt' ,fignumber )
fignumber=fignumber+1;
%%
disp('Away from Global Minimum and surrounded by  stationary points around')
x0=[0,-4]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=steepest_optim_ls(@func_basis,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', xopt' ,fignumber )
fignumber=fignumber+1;
%%
disp('Close from Global Minimum and surrounded by  stationary points around')
x0=[3,-0.5]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=steepest_optim_ls(@func_basis,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', xopt' ,fignumber )
fignumber=fignumber+1;
%% Q2.4:  Newton's algorithm

fignumber=1;

x_global_min=[3,2];

disp('Pretty Close to the global minimum and far away from any other stationary point')
x0=[5,5]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=newton_optim_ls(@func_basis,@func_basis_2,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', xopt' ,fignumber )
fignumber=fignumber+1;
%%
disp('Close to Global Minimum and surrounded by stationary points')
x0=[0,0]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=newton_optim_ls(@func_basis,@func_basis_2,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', xopt' ,fignumber )
fignumber=fignumber+1;
%%
disp('Far away from Global Minimum and surrounded by stationary points')
x0=[-5,0]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=newton_optim_ls(@func_basis,@func_basis_2,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', xopt' ,fignumber )
fignumber=fignumber+1;
%%
disp('Away from Global Minimum and surrounded by  stationary points around')
x0=[0,-4]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=newton_optim_ls(@func_basis,@func_basis_2,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', xopt' ,fignumber )
fignumber=fignumber+1;
%%


disp(' Very Close from Global Minimum and surrounded by  stationary points around')
x0=[3,-0.5]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=newton_optim_ls(@func_basis,@func_basis_2,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', xopt' ,fignumber )
fignumber=fignumber+1;
%% Q2.5:BFGS



fignumber=1;

x_global_min=[3,2];

disp('Pretty Close to the global minimum and far away from any other stationary point')
x0=[5,5]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=bfgs_optim_ls(@func_basis,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', xopt' ,fignumber )
fignumber=fignumber+1;
%%
disp('Close to Global Minimum and surrounded by stationary points')
x0=[0,0]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=bfgs_optim_ls(@func_basis,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', xopt' ,fignumber )
fignumber=fignumber+1;
%%
disp('Far away from Global Minimum and surrounded by stationary points')
x0=[-5,0]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=bfgs_optim_ls(@func_basis,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', xopt' ,fignumber )
fignumber=fignumber+1;
%%
disp('Away from Global Minimum and surrounded by  stationary points around')
x0=[0,-4]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=bfgs_optim_ls(@func_basis,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', xopt' ,fignumber )
fignumber=fignumber+1;
%%
disp('Close from Global Minimum and surrounded by  stationary points around')
x0=[3,-0.5]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=bfgs_optim_ls(@func_basis,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', xopt' ,fignumber )
fignumber=fignumber+1;
