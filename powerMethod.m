function [eigvals, eigvecs] = powerMethod(A)
    % Computes all eigenvalues and eigenvectors using
    %
    n = size(A,1);
    eigvals = zeros(n,1);
    eigvecs = zeros(n,n);
    B = A; 

    for k = 1:n
        x = rand(n,1);
        x = x / norm(x);
        for iter = 1:500
            y = B * x;
            x_new = y / norm(y);
            lambda = (x_new' * B * x_new) / (x_new' * x_new);
            if norm(x_new - x) < 1e-12
                x = x_new;
                break;
            end
            x = x_new;
        end
        eigvals(k) = lambda;
        eigvecs(:,k) = x;
        B = B - lambda * (x * x');
    end
end


