function [x_star,r_star] = NOfit_Phil(t,y,n)
%Saving the y-data from assignment here for easy sake
Y = [91.2200 28.0400 22.9100 26.6500 42.9600 101.0500 202.3600 328.0200 364.1200 299.2300 238.0000 227.4900 218.0300 223.6200 238.7500 271.2600 267.7200 251.3200 230.0400 206.6900 170.7700 131.6700 143.8500 157.5700];

%Converts data vector y into vector (Nx1) format
if size(y,1) == 1
    y = y';
end

if mod(n,2) ~= 1
    x_star = 'You must pick an odd value for the fit order';
    r_star = [];

else
    % Initial frequency constant
    w = 2*pi/24;

    % Generating function value matrix A
    A= zeros(length(t),n);
    A(:,1) = 1;
        
        for i = 2:2:(n-1)
            N = i/2;
            A(:,i) = sin(w*N*t);
            A(:,i+1) = cos(w*N*t);
        end
        
        x_star = A\y;
        r_star = y - A*x_star;
        
        % Plotting the data points against the model
        
        %Creating a more extensive function value matrix B for a smooth
        %model plot
        
        t_new = 1:0.1:24;
        B = zeros(length(t_new),n);
        B(:,1) = 1;
        
        for i = 2:2:(n-1)
            N = i/2;
            B(:,i) = sin(w*N*t_new);
            B(:,i+1) = cos(w*N*t_new);
        end
        
        figure(1)
        s = scatter(t,y,50,'filled');
        s.MarkerFaceColor = 'r';
        hold on
        p = plot(t_new,B*x_star,'-');
        p.LineWidth = 1.2;
        p.Color = 'b';
        grid on
            
end
end

