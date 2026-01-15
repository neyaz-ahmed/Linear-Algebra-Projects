% Gaussian Elimination Method
% No of equation = No of unknown
%%
function x = gauss(A,b)

Ab = [A b];       % Augmented matrix
n = size(A,1);

% Forward Elimination
for k = 1:n-1
    % Pivoting (if needed, for stability)
    [~, idx] = max(abs(Ab(k:n,k)));
    idx = idx + k - 1;
    
    if idx ~= k                % Swap rows        
      Ab([k idx], :) = Ab([idx k], :);
    end
    
    for i = k+1:n
        factor = Ab(i,k)/Ab(k,k);
        Ab(i,k:end) = Ab(i,k:end) - factor*Ab(k,k:end);
    end

end

% Back Substitution to find unknowns
x = zeros(n,1);
x(n) = Ab(n,end)/Ab(n,n);

for i = n-1:-1:1
    x(i) = (Ab(i,end) - Ab(i,i+1:n)*x(i+1:n)) / Ab(i,i);
end


end