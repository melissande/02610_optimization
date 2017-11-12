function  convergence_proof( df, x, xopt ,fignumber )
%convergence proof plots the infinite norm of the gradient of f func
%and the deviation of the estimation x from the global minimum xopt

%% Arguments:
%df:gradient vector  taken at each estimation point for each iteration (size: number of iterations x 2)
%x: vector of estimation at each iteration(size: number of iterations x 2)
%xopt: minimizer estimtated at the last iteration (size: 1x 2)
%fignumber: dummy (considered for display )
iter=(1:length(x))';
figure(fignumber)
subplot(2,1,1)
plot(iter,max(abs(df),[],2),'ok')
xlabel('Iterations')
ylabel('||df_k||_{inf}')
subplot(2,1,2)
plot(iter,sqrt(sum((x-xopt).^2,2)),'ok')
xlabel('Iterations')
ylabel('e_k')
end

