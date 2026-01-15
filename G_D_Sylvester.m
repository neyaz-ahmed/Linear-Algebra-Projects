%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sylvester equation: A*X + X*B = C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all;

%% 
n=20;
A=rand(n,n);    C=eye(n);   B=eye(n);
%%
n = size(A,1);
K = kron(eye(n), A) + kron(B', eye(n));
b = reshape(C, n^2, 1);
vecX = K \ b;
X = reshape(vecX, n, n);

X_exact = sylvester(A, B, C);
fprintf('||X - Xexact||_F = %g (should be ~0)\n', norm(X - X_exact, 'fro'));


