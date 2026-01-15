%% MAtlab SCript for Quadratic polynomial fitting
clear all;
clc;

%% Load your data here
x = [-1, 0, 1, 2, 3];
y = [-6.2; -2; 0; 1.5; -2]';
x=x + rand(1,5);
%%  Construct system of linear equations:  AX=y
A = [x.^2; x; ones(1, length(x))]';

%%  Solve  AX=y
X = (A'*A) \ (A'*y');

% Evaluate the polynomial at the given x points
y_fit = X(1)*x.^2 + X(2)*x + X(3);

%% Plot the approximate polynomial
xt = linspace(min(x)-1, max(x)+1, 100);
y_approx = X(1)*xt.^2 + X(2)*xt + X(3);

%% Plot the given points
figure(1);
plot(x, y, 'ro', 'LineWidth', 1.5); 
hold on;
plot(xt, y_approx, 'b-.', 'LineWidth', 2); 
xlabel('x');
ylabel('y');
title('Example 1: Quadratic Polynomial Fit');
legend('Given Points', 'Fitted Quadratic');
grid on;
hold off;

%% display the sum of squared residuals
residuals = y' - y_fit;
err = norm(residuals, 2);
fprintf('Error norm: %.4f\n', err);

