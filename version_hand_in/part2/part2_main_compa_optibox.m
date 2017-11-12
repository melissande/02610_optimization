%% PART 2  - Compa OPTIBOX

%% Function definition

func=@(x1,x2)((x1.^2+x2-11).^2+(x1+x2.^2-7).^2);

%% Q2.1: Contour Plots

%Definition of stationary points (found analytically)
x_sad=[ 0.08667,2.88425;-3.07302,-0.08135; 3.38515,0.07385;-0.12796, -1.95371];
x_min=[ -2.80511,3.13131;3.58442,-1.84812;-3.77931,-3.28318];
x_max=[-0.27084,-0.92303];
x_global_min=[3,2];


%% ucfmin
fignumber=1;

x_global_min=[3,2];

disp('Pretty Close to the global minimum and far away from any other stationary point')
x0=[5,5]';
fprintf('x0=[%d,%d]\n',x0);
[X, info,perf,tmp] = ucminf(@func_basis, x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated \n',X(:,end),info(4),info(5));
error=sqrt(sum((X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( perf.ng', X', x_global_min ,fignumber )
fignumber=fignumber+1;
%%
disp('Close to Global Minimum and surrounded by stationary points')
x0=[0,0]';
fprintf('x0=[%d,%d]\n',x0);
[X, info,perf,tmp] = ucminf(@func_basis, x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated \n',X(:,end),info(4),info(5));
error=sqrt(sum((X'-x_global_min).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( perf.ng', X', x_global_min ,fignumber )
fignumber=fignumber+1;
%%
disp('Far away from Global Minimum and surrounded by stationary points')
x0=[-5,0]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=bfgs_my_ls(@func_basis,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_min(1,:)).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', x_min(1,:) ,fignumber )
fignumber=fignumber+1;
%%
disp('Away from Global Minimum and surrounded by  stationary points around')
x0=[0,-4]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=bfgs_my_ls(@func_basis,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_min(3,:)).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min,fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', x_min(3,:) ,fignumber )
fignumber=fignumber+1;
%%
disp('Close from Global Minimum and surrounded by  stationary points around')
x0=[3,-0.5]';
fprintf('x0=[%d,%d]\n',x0);
[xopt,stat]=bfgs_my_ls(@func_basis,x0);

fprintf('Global Minimum estimated: [%d,%d] in %d iterations,%d functions evaluated in %d seconds\n',xopt,stat.iter,stat.nfun,stat.tmp);
error=sqrt(sum((stat.X'-x_min(2,:)).^2,2));
[p,c]= cvg_rate( error, fignumber);
fignumber=fignumber+1;
fprintf('Convergence rate: %d and constant limit : %d\n',p,c);
make_contour_plot( func, stat.X', x_global_min',fignumber,[] )
fignumber=fignumber+1;
convergence_proof( stat.dF', stat.X', x_min(2,:) ,fignumber )
fignumber=fignumber+1;

