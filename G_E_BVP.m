
clear all;
clc;

[A,b,x] = data_BVP();

%% Solve using LU decomposition
[L, U] = lu(A);
% y = L \ (P * b); 
% u_lu = U \ y;   
u_lu = L * U \ b;

disp('LU: Residual norm ||Au - b||_2:');
disp(norm(A * u_lu - b));

%% Solve using QR decomposition
[Q, R] = qr(A);
% u_qr = R \ (Q' * b); 
u_qr = Q * R \ b;
disp('QR: Residual norm ||Au - b||_2:');
disp(norm(A * u_qr - b));
%% SVD 

[U, S, V] = svd(A);

u_svd = U * S * V' \ b;
disp('SVD: Residual norm ||Au - b||_2:');
disp(norm(A * u_svd - b));
%% Gauss Elimination 
u_g = gauss(A,b); 

disp('Gauss: Residual norm ||Au - b||_2:');
disp(norm(A * u_g - b));
%% %% Jacobi Iteration
X_jac = jacobi(A, b);

%% Gauss Seidel
X_gs = gauss_seidel(A, b);

disp('Gauss Seidel: Residual norm ||Au - b||_2:');
disp(norm(A * X_gs - b));
%% Plot the solution
figure(11);
plot(x, u_lu, 'b-', 'LineWidth', 1, 'DisplayName', 'LU Solution');
hold on;
plot(x, u_qr, 'r-.', 'LineWidth', 1, 'DisplayName', 'QR Solution');
plot(x, u_svd, 'r-.', 'LineWidth', 1, 'DisplayName', 'SVD Solution');

xlabel('x');
ylabel('u(x)');
title('Finite Difference Solution to Boundary Value Problem');
legend('show');
grid on;

