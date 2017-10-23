function [] = make_contour_plot_init( func, x_sad,x_min,x_max, xopt )
%make_contour_plot build a contour plot of func which should be
% defined as an anonymous function and plots points where the gradient is 0
%(saddle points, local minima, local maxima) which had been found
%anatically
%% INPUT
%func: anonymous function
%x_sad:saddle points coordinates
%x_min:local minima coordinates
%x_max:local maxima coordinates
%x_opt:global minimum
%% OUTPUT
%figure with the contour plot and theoretical and algorithms minimum


x1=linspace(-5,5,1000);
x2=linspace(-5,5,1000);

[X1,X2]=meshgrid(x1,x2); 
func_eval=func(X1,X2);

% Make contour plot
figure;
fs = 10; 
contour(X1,X2,func_eval,100,'linewidth',1)

iter_sad=1:size(x_sad,1);
if ~isempty(x_sad)
    hold on
    scatter(x_sad(:,1), x_sad(:,2),30,'black','filled')
    b = num2str(iter_sad'); c = strcat('saddle ',cellstr(b));
    dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
    text(x_sad(:,1)+dx, x_sad(:,2)+dy, c,'Color','black','FontSize',14);
end
iter_min=1:size(x_min,1);
if ~isempty(x_min)
    hold on
    scatter(x_min(:,1), x_min(:,2),30,'black','filled')
    b = num2str(iter_min'); c = strcat('min ',cellstr(b));
    dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
    text(x_min(:,1)+dx, x_min(:,2)+dy, c,'Color','black','FontSize',14);
end

iter_max=1:size(x_max,1);
if ~isempty(x_max)
    hold on
    scatter(x_max(:,1), x_max(:,2),30,'black','filled')
    b = num2str(iter_max'); c = strcat('max ',cellstr(b));
    dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
    text(x_max(:,1)+dx, x_max(:,2)+dy, c,'Color','black','FontSize',14);
end

if ~isempty(xopt)
    hold on
    scatter(xopt(1),xopt(2),100,'blue','filled','h')
    dx = 0.1; dy =0.1; % displacement so the text does not overlay the data points
    text(xopt(1)+dx, xopt(2)+dy, 'Min','Color','blue','FontSize',14);
end
xlabel('x_1','fontsize',fs); 
ylabel('x_2','fontsize',fs); 
set(gca,'fontsize',fs);
axis equal
grid on
title('Contour Plot')
colormap hsv;
colorbar

end