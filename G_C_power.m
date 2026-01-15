%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute eigenvalues and eigenvectors of a matrix
clear all;
clc; 
% n = 5;             
% A = triu( rand( n,n ) );
%%
n = 200; m = 10; epsilon = 1e-12;
A = spdiags([-ones(n,1) 2*ones(n,1) -ones(n,1)], -1:1, n, n);

 %%  power iteration
[eig2, vec2] = powerMethod(A);    
disp('Power: Check A*v -v*lambda:');
disp(norm(A * vec2 - vec2 * diag(eig2)));
 %%   Inverse iteration with Rayleigh quotient shifts.
[eig3, vec3, iter] = inversepower(A);    
disp('Inverse: Check A*v -v*lambda:');
disp(norm(A * vec3 - vec3 * diag(eig3)));
 
%%   use matlab eig(A)
[V,D]=eig(full(A));
fprintf('Original eigen values: eig(A) \n \n');
eig1=sort(diag(D));
eig2=sort(eig2);
eig3=sort(eig3);
disp(norm(eig1-eig2))
disp(norm(eig1-eig3))

