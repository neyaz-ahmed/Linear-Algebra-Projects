% Electrical Circuit Analysis: RLC Circuit Solution
clear; clc; 
close all;
%% Load data
R = 1;      
L = 1;     
C = 0.5;    

A = [0,    1/L; 
    -1/C, -1/(R*C)];

[V, D] = eig(A);
lambda = diag(D);

%% 
fprintf('\n=== PART (b): Specific Solution with Initial Conditions ===\n');
% Initial conditions
iL0 = 2;    
vC0 = 1;    
x0 = [2; 1];
c = V \ x0;
fprintf('Constants: c1 = %.4f, c2 = %.4f\n', c(1), c(2));

%% Numerical Solution 
fprintf('\n=== Numerical Solution and Visualization ===\n');
t = linspace(0, 3, 50);
iL_t = c(1)*V(1,1)*exp(lambda(1)*t) + c(2)*V(1,2)*exp(lambda(2)*t);
vC_t = c(1)*V(2,1)*exp(lambda(1)*t) + c(2)*V(2,2)*exp(lambda(2)*t);

%% Plot results
figure(10);
subplot(2,1,1);
plot(t, iL_t, 'b-', 'LineWidth', 2);
grid on;
title('Inductor Current i_L(t)');
xlabel('Time t (seconds)');
ylabel('i_L(t) (amperes)');
legend('i_L(t)');

subplot(2,1,2);
plot(t, vC_t, 'r-', 'LineWidth', 2);
grid on;
title('Capacitor Voltage v_C(t)');
xlabel('Time t (seconds)');
ylabel('v_C(t) (volts)');
legend('v_C(t)');

