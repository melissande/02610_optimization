function  convergence_proof( df, x, xopt ,fignumber )
%convergence proof plots the infinite norm of the gradient of f func
%and the deviation of the estimation x from the global minimum xopt



iter=(1:length(x))';
figure(fignumber)
subplot(2,1,1)
plot(iter,max(df,[],2),'ok')
xlabel('Iterations')
ylabel('||df_k||_{inf}')
subplot(2,1,2)
plot(iter,sqrt(sum((x-xopt).^2,2)),'ok')
xlabel('Iterations')
ylabel('e_k')
end

