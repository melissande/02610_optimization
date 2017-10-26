function [x,stat] = GaussNewton(fun,x0)

% Solver settings and info
tic;

maxit = 100*length(x0);
tol   = 1.0e-10;

stat.converged = false; % converged
stat.nfun      = 0;     % number of function calls
stat.iter      = 0;     % number of iterations
stat.tmp       = 0;     %time stamp
stat.nJac      = 0;     % number of jacobian calls

% Initial iteration
x = x0;
it = 0;

%Initial function computation
[f,df,J] = feval(fun,x);

converged = (norm(df,'inf') <= tol);
stat.nfun = 1;
stat.nJac = 1;

% Store data for plotting
stat.X  = x;
stat.F  = f;
stat.dF = df;

% Main Loop of Gauss Newton
while ~converged && (it < maxit)
    it = it+1;
    
    % Gauss Newton Algorithm with Line Search
    % ================================================
    
    d = - (J'*J)\(df);

    [x,~,~,eval] = my_line_search(fun,x,f,df,d,stat.nfun,[]);
    
    stat.nfun = eval;
    
    % ================================================
    
    [f,df,J] = feval(fun,x);
    
    converged = (norm(df,'inf') <= tol);
    
    stat.nfun = stat.nfun +1;
    stat.nJac = stat.nJac +1;
    
    % Store data for plotting

    stat.X  = [stat.X  x];
    stat.F  = [stat.F f];
    stat.dF = [stat.dF df];
end

% Prepare return data
if ~converged
x = []; 
end

stat.converged = converged;
stat.iter = it;

toc;
stat.tmp = toc;
end
