%% Large-scale Balanced Truncation Example (corrected)
clc; 
clear all;
close all;

%% Load System matrices
n = 20; m = 2; p = 1;
rng(1);
A = randn(n); 
A = -A*A' - 0.5*eye(n); 
B = ones(n, m);
C = 2*ones(p, n);

%% Compute Controllability and Observability Gramians
Wc = lyap(A, B*B'); 
Wo = lyap(A', C'*C);

% Compute Hankel Singular Values (Singular Values of sqrt(Wc)*sqrt(Wo))
[U, S, V] = svd(Wc * Wo);
r=5;
% Compute Transformation Matrices
U_r = U(:, 1:r);  
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);
% Balanced transformation
T = (S_r)^(-1/2) * V_r' * Wo;
T_inv = Wc*U_r*(S_r)^(-1/2);

% Reduced System Matrices

Ar =T* A * T_inv;
Br = T*B ;
Cr = C * T_inv;
% D_r = D;  % Remains unchanged

%% Step 6: Compare step responses
sys_full = ss(A,B,C,0);
sys_red = ss(Ar,Br,Cr,0);
t = 0:0.1:2;
u = ones(length(t),1); % Unit step input for channel 1
y_full = lsim(sys_full(1,1), u, t);
y_red = lsim(sys_red(1,1), u, t);

% Step 7: Print Hankel singular values
hsv = diag(S);
% Step 7.1: Hankel singular value comparison (plot)
figure(4);
semilogy(1:length(hsv), hsv, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
semilogy(r, hsv(r), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Index'); ylabel('Hankel Singular Value');
title('Hankel Singular Values (log scale)');
legend('HSVs', sprintf('Truncation at r=%d', r), 'Location', 'best');
grid on;
% Step 8: Error estimate
error_bound = 2*sum(hsv(r+1:end));
fprintf('Error bound (2*sum of neglected HSVs) = %g\n', error_bound);

% Step 9: Relative error comparison
rel_error_l2 = norm(y_full - y_red) / norm(y_full);
rel_error_max = max(abs(y_full - y_red) ./ abs(y_full + eps)); % Avoid div by zero
fprintf('Relative L2 error: %g\n', rel_error_l2);
fprintf('Maximum relative error: %g\n', rel_error_max);

% Step 10: Plot relative error over time
rel_error_time = abs(y_full - y_red) ./ abs(y_full + eps);
figure(5);
plot(t, rel_error_time, 'g-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Relative Error');
title('Pointwise Relative Error Over Time');
grid on;

figure(3);
plot(t, y_full, 'b', 'LineWidth', 1.5); hold on;
plot(t, y_red, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Output');
% title('Balanced Truncation: Original vs Reduced');
legend('Original','Reduced');
grid on;