% Implement Arnoldi iteration for a 4x4 matrix
clear all;
clc;
%%
n = 10; m = n;
epsilon = 1e-8;
% A = 0.01*randi(n,n); % 
A=spdiags([-ones(n,1) 2*ones(n,1) -ones(n,1)], -1:1, n, n);
A=A'*A;
b = randn(n,1);    b = b/norm(b); 
%%  Arnoldi iteration 
tic; [V, H] = arnoldi(A, b, m); toc;
fprintf('Arnoldi: Orthonormality check (V^T V):\n'); disp(norm(V' * V))

fprintf('Verification: norm(AV - V H):\n')
disp(norm(A * V(:,1:m) - V * H))

%% Lanczos iteration (symmetric case)
tic;[V, T] = lanczos(A, b, m);toc;
fprintf('Lanczos: Orthonormality check (V^T V):\n'); disp(norm(V' * V))

fprintf('Verification: norm(AV - V H):\n')
disp(norm(A * V(:,1:m) - V * T))
%%  Householder reflection 
 [V, H] = householder(A); 
 fprintf('Lanczos: Orthonormality check (V^T V):\n'); disp(norm(V' * V))

fprintf('Verification: norm(AV - V H):\n')
disp(norm(A * V(:,1:m) - V * H))

 
 
