function [x,stat] = bfgs_my_ls(fun,x0)
% Solver settings and info
tic;
%maxit=50;
maxit = 100*length(x0);
tol   = 1.0e-5;
B  = eye(length(x0));
stat.converged = false;% converged
stat.nfun      = 0;% number of function calls
stat.iter      = 0;% number of iterations


% Initial iteration
x = x0;
it = 0;



[f,df] = feval(fun,x);

converged = (norm(df,'inf') <= tol);
stat.nfun = 1;
stat.tmp=0;
% Store data for plotting
stat.X = x;
stat.F = f;
stat.dF = df;

% Main loop of steepest descent
while ~converged && (it < maxit)
    it = it+1;
    % Newton's algo
    % TODO -- Insert code between the lines
    % ================================================
    d = -B \ df;

    [xnew,~,~,eval] =my_line_search(fun, x,f,df, d,stat.nfun,[]);
    
    stat.nfun=eval;
    
    
    % ================================================
    [fnew,dfnew] = feval(fun,xnew);
    
    
    converged = (norm(dfnew,'inf') <= tol);
    
    s = xnew - x;
    y = dfnew - df;
    Bs = B*s;
    B = B - (Bs/(s'*Bs))*Bs'  + (y/(y'*s))*y';
    x = xnew;
    f = fnew;
    df = dfnew;
    stat.nfun = stat.nfun+1;
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
stat.tmp=toc;
end