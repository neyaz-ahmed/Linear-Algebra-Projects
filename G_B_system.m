% Solving a Hilbert system using different decompositions
% -------------------------------------------------------------
clear all; clc;

%% Load your data for AX=B
% n = 20;
% A = rand(n);
% A = floor(10*rand(n)); 
%%
n = 50;
A = zeros(n, n);
A= rand(n,n);
for i = 1:n
    A(i,i) = 400; 
    if i > 1
        A(i,i-1) = -1; 
        A(i-1,i) = -1; 
    end
end
x_exact=floor(10*rand(n,1));              

%% Choose exact solution 
% x_exact = ones(n,1);      
%x_exact = rand(n,1); 
% spy(A)
B = A * x_exact;           % Construct the right-hand side
% noise = rand(n,1);
% B= B + noise;
%%  Gauss elimination 
tic; X_g = gauss(A,B);   toc; 

%% LU decomposition 
tic; X_lu = solveLU(A,B); toc;
 
%% QR decomposition   
tic;  [Q, R] = qr(A, 0); 
X_qr = R \ (Q' * B);  toc; 

%% SVD  decomposition  
tic;  [U, S, V] = svd(A, 'econ'); 
X_svd = V * (S \ (U' * B));  toc; 

%%  Cholesky decomposition
Aspd=A'*A; b=A'*B;
tic; X_ch = solvechol(Aspd,b); toc;

%% Jacobi Iteration
X_jac = jacobi(A, B);

%%   Gauss-Seidal 
X_gs = gauss_seidel(A, B);
%%  CG
A1=A'*A; B1=A'*B;  x0 = zeros(n, 1); tol = 1e-10; 
X_cg = conjugate_gradient(A1, B1, x0, tol);

%%    GMRES
X_gm = gmres(A,B,x0);


%% exact solution
tic; X_b = A \ B; toc;

%%  Errors
fprintf('Errors for Gaussian elimination:    %.2e\n', norm(x_exact - X_g));
fprintf('Errors for LU decomposition:        %.2e\n', norm(x_exact - X_lu));
fprintf('Errors for QR decomposition:        %.2e\n', norm(x_exact - X_qr));
fprintf('Errors for SVD:                     %.2e\n', norm(x_exact - X_svd));
fprintf('Errors for Cholesky decomposition:  %.2e\n', norm(x_exact - X_ch));
fprintf('Errors for Jacobi:                  %.2e\n', norm(x_exact - X_jac));
fprintf('Errors for Gauss-Seidal:            %.2e\n', norm(x_exact - X_gs));
fprintf('Errors for GMRES:                   %.2e\n', norm(x_exact - X_gm));
fprintf('Errors for CG:                      %.2e\n', norm(x_exact - X_cg));
%% condition number using infinity norm
cond_inf = cond(A, inf);
fprintf('\n Condition number of A is %.2e \n', cond_inf);


