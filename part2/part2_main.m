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
    [f,df,d2f]=func_basis_2(x_sad(i,:));
   fprintf('Saddle point %d: df=[%d,%d]  and d2f=[%d,%d;%d,%d] \n',i,df,d2f)
end

for i=1:size(x_min,1)
    [f,df,d2f]=func_basis_2(x_min(i,:));
   fprintf('Local minimum %d: df=[%d,%d]  and d2f=[%d,%d;%d,%d] \n',i,df,d2f)
end

for i=1:size(x_max,1)
    [f,df,d2f]=func_basis_2(x_max(i,:));
   fprintf('Local maximum %d: df=[%d,%d]  and d2f=[%d,%d;%d,%d] \n',i,df,d2f)
end

[f,df,d2f]=func_basis_2(x_global_min);
fprintf('Global minimum %d: df=[%d,%d]  and d2f=[%d,%d;%d,%d] \n',i,df,d2f)


%% Q2.4:  Newton's algorithm
x0=[3.1,2.1]';
[xopt,stat]=newton_my_ls(@func_basis,@func_basis_2,x0);

