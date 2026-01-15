%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Simple MATLAB script for computing e^At 
clear all; clc;
%% Define the matrix A
A = [-3 1; 1 -4];
[V,D]=eig(A);
 S =V; 
 S_inv =inv(S); 

syms t;  
eDt = [exp(D(1)*t)     0; 
       0           exp(D(2)*t)];
eAt = S * (eDt * S_inv);

f = matlabFunction(eAt);

t_val = 1;
eAt_numeric_our = f(t_val);

eAt_numeric_matlab = expm(A * t_val);
fprintf('MATLAB built-in expm(A*%d):\n', t_val);
disp(eAt_numeric_matlab);

difference = norm(eAt_numeric_our - eAt_numeric_matlab, 'fro');
fprintf('Difference between methods: %e\n', difference);

