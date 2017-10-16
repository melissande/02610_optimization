function [] = make_contour_plot( func, x, xopt )
%make_contour_plot build a contour plot of func which should be
% defined as an anonymous function and plots the minimum at each step (x)
% of the algorithm and the theoretical minimum (xopt)
%% INPUT
%func: anonymous function
%x: vector (or matrix if dim>1)with all the minimum found at each step of the algorithm (can be empty)
%xopt: scalar (or vector if dim>1)  with the theoretical minimum
%% OUTPUT
%figure with the contour plot and theoretical and algorithms minimum
iter=1:length(x);

x1=linspace(-100,100,100);
x2=linspace(-100,100,100);

[X1,X2]=meshgrid(x1,x2); 
func_eval=func(X1,X2);

% Make contour plot
figure;
fs = 10; 
contour(X1,X2,func_eval,50,'linewidth',2)
if ~isempty(iter)
    hold on
    scatter(x(1,:), x(2,:),30,'red','filled')
    b = num2str(iter'); c = cellstr(b);
    dx = 2; dy = 2; % displacement so the text does not overlay the data points
    text(x(1,:)+dx, x(2,:)+dy, c,'Color','red','FontSize',14);
end
if ~isempty(xopt)
    hold on
    scatter(xopt(1),xopt(2),100,'green','filled','h')
    dx = 2; dy =2; % displacement so the text does not overlay the data points
    text(xopt(1)+dx, xopt(2)+dy, 'Min','Color','green','FontSize',14);
end
xlabel('x_1','fontsize',fs); 
ylabel('x_2','fontsize',fs); 
set(gca,'fontsize',fs);
axis equal
grid on
title('Contour Plot')
colorbar

end

