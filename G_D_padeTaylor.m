%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
clc;
%%   Load data
n = 10;                 
A = triu( rand( n,n ) );
I = eye(size(A, 1));

%%   Taylor series (order 3)
A2 = A * A;
A3 = A2 * A;
taylor_A = I + A + A2/2 + A3/6;

%%   Pade approximation (2,2)
P = I + 0.5 * A + A2/12;
Q = I - 0.5 * A + A2/12;
pade_A = Q \ P;

%%  Using eigen decomposition
[V, D] = eig(A);           
fD = diag(exp(diag(D)));     
fA_eig = V * fD * inv(V);        

%%    Using Schur decomposition
[Q, U] = schur(A, 'real');  
fT = expm(U);                
fA_schur = Q * fT * Q';     

%%     Jordan decomposition
[V, J] = jordan(A);                
fJ = zeros(size(J));
for k = 1:size(J,1)
    if k > 1 && J(k,k-1) ~= 0       
        fJ(k,k-1) = exp(J(k,k));    
    else
        fJ(k,k) = exp(J(k,k));
    end
end
fA_jordan = V * fJ * inv(V); 

%% Krylov subspace method (for e^A * b)
m = 3; 
b = rand(size(A, 1),1);
b=b/ norm(b); 
[V, H] = arnoldi(A, b, m);
krylov_A_b = V(:,1:m) * expm(H(1:m,1:m)) * (V(:,1:m)' * b);

%% Exact computation by MATLAB
exact_A = expm(A);
%%    Error testing
error_taylor = norm(exact_A - taylor_A, 'fro');
error_pade = norm(exact_A - pade_A, 'fro');
error_eig = norm(exact_A - fA_eig, 'fro');
error_schur = norm(exact_A - fA_schur, 'fro');
error_jordan = norm(exact_A - fA_jordan, 'fro');
exact_A_b = exact_A * b;
error_krylov = norm(exact_A_b - krylov_A_b, 2);

fprintf('Error Taylor: %.4e\n', error_taylor);
fprintf('Error Pade: %.4e\n', error_pade);
fprintf('Error eig: %.4e\n', error_eig);
fprintf('Error schur: %.4e\n', error_schur);
fprintf('Error Jordan: %.4e\n', error_jordan);
fprintf('Error Krylov (e^A*b): %.4e\n', error_krylov);


