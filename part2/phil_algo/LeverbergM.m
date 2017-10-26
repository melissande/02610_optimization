function [x,stat] = LeverbergM(fun,x0)

% Solver settings and info
tic;

%Stopping Constraints
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
tau = 1.0e-3;

%Initial function computation
[f,df,J,r] = feval(fun,x);

converged = (norm(df,'inf') <= tol);
stat.nfun = 1;
stat.nJac = 1;

%Initial uI computation
A = J'*J;
I = eye(length(diag(A)));
u = tau*max(diag(A));

% Store data for plotting
stat.X  = x;
stat.F  = f;
stat.dF = df;

% Main Loop of Gauss Newton
while ~converged && (it < maxit)
    it = it+1;
    
    % Leverberg-Marquandt Algorithm
    % ================================================
    
    d = -(J'*J+u*I)\df;
    
    % ================================================
    
    %Additional Stopping Constraint
    if norm(d) <= tol*(norm(x)+tol)
        converged = true;  
    else
        %Update uI through p calculation
        x_new = x+d;
        [f,df_new,J,r_new] = feval(fun,x_new);
        
        stat.nfun = stat.nfun +1;
        stat.nJac = stat.nJac +1;
        
        p = ((1/2)*(r-r_new)'*(r-r_new))/((1/2)*d'*(u*d-df));
        
        if p > 0
            u = u*max([1/3,1-(2*p-1)^3]);
            converged = (norm(df_new,'inf') <= tol);
        else
            u = 2*u;
        end
        
        %Update x, df and r
        x  = x_new;
        df = df_new;
        r  = r_new;
    
        % Store data for plotting

        stat.X  = [stat.X  x];
        stat.F  = [stat.F f];
        stat.dF = [stat.dF df];
    end
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
