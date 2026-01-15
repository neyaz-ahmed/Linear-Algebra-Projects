%% Solve Lyapunov using 9x9 linear system and backslash
clc; clear all;
  
n = 20; 
rng(1);
A = floor(randn(n)) - 2*eye(n); 
B = ones(n, n);
C = 2*ones(n, n);

K = kron(eye(n), A) + kron(A, eye(n));   
b = -reshape(C, n^2, 1); 
vecX = K \ b;
X = reshape(vecX, n, n);

disp('Solution X (from 9x9 solve):');
disp(X);

residual = norm(A*X + X*A' + C);   % should be (close to) zero
disp('Residual A*X + X*A'' + C:');
disp(residual);

X_exact = lyap(A, C);
fprintf('||X - Xexact||_F = %g  (should be ~0)\n', norm(X - X_exact, 'fro'));

