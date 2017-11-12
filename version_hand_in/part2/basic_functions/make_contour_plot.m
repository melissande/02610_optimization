function [] = make_contour_plot( func, x, xopt ,fignumber,zoom)
%make_contour_plot build a contour plot of func which should be
% defined as an anonymous function and plots the minimum at each step (x)
% of the algorithm and the theoretical minimum (xopt)
%% INPUT
%func: anonymous function
%x: vector (or matrix if dim>1)with all the minimum found at each step of the algorithm (can be empty)
%xopt: scalar (or vector if dim>1)  with the theoretical minimum
%% OUTPUT
%figure with the contour plot and theoretical and algorithms minimum
iter=1:size(x,1);

x1=linspace(-5,5,500);
x2=linspace(-5,5,500);

[X1,X2]=meshgrid(x1,x2); 
func_eval=func(X1,X2);

% Make contour plot
figure(fignumber);
fs = 10; 
contour(X1,X2,func_eval,50,'linewidth',1)
if ~isempty(iter)
    hold on
    plot(x(:,1), x(:,2),'-ok')
    %scatter(x(:,1), x(:,2),30,'black','filled')
    %b = num2str(iter'); c = cellstr(b);
    %dx = 0.05; dy =0.05;  % displacement so the text does not overlay the data points
    %text(x(:,1)+dx, x(:,2)+dy, c,'Color','black','FontSize',14);
    
end
if ~isempty(xopt)
    hold on
    scatter(xopt(1),xopt(2),100,'green','filled','h')
    
end
xlabel('x_1','fontsize',fs); 
ylabel('x_2','fontsize',fs); 
set(gca,'fontsize',fs);
grid on
if ~isempty(zoom)
    xlim(zoom(1,:))
    ylim(zoom(2,:))
else 
    axis equal
colormap hsv;
colorbar

end

