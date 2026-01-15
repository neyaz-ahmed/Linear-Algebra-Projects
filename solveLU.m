function x = solveLU(A,b)

% MATLAB's lu function returns L, U, and permutation matrix P such that P*A = L*U
% Ax=b    => PAx=Pb   => LUx=Pb   => x=U\(L\(Pb)) 

[L, U, P] = lu(A); 


y = L \ (P * b);

x = U \ y;

% disp('Verification: Norm of P*A - L*U (should be close to 0):');
% disp(norm(P*A - L*U));

end