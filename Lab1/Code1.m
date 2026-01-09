clc;        % Clear command window
clear;      % Clear workspace
close all;  % Close all figures

% -----------------------------
% Given data
% -----------------------------
t = 1:10;   % t_i values
y = [1.3 3.5 4.2 5.0 7.0 8.8 10.1 12.5 13.0 15.6];  % y_i values

% -----------------------------
% Fine grid for smooth plotting
% -----------------------------
t_fine = linspace(1,10,200);

% -----------------------------
% Create a figure
% -----------------------------
figure;
hold on;
grid on;

% Plot original data
plot(t, y, 'ro', 'MarkerSize', 8, 'LineWidth', 2);

% -----------------------------
% Loop for degrees 1 to 4
% -----------------------------
for n = 1:4
    
    % Least squares polynomial fit
    p = polyfit(t, y, n);
    
    % Evaluate polynomial at data points
    y_fit = polyval(p, t);
    
    % Compute error E (sum of squared errors)
    E = sum((y - y_fit).^2);
    
    % Display results
    fprintf('\nDegree %d polynomial coefficients:\n', n);
    disp(p);
    fprintf('Error E = %.6f\n', E);
    
    % Plot polynomial curve
    plot(t_fine, polyval(p, t_fine), 'LineWidth', 2);
end

% -----------------------------
% Labels and legend
% -----------------------------
xlabel('t');
ylabel('y');
title('Least Squares Polynomial Fits (Degree 1 to 4)');
legend('Data', 'Degree 1', 'Degree 2', 'Degree 3', 'Degree 4', 'Location', 'northwest');
