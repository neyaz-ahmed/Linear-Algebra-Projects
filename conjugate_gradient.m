 function x = conjugate_gradient(A, b, x0, tol)
% Reminder: A should be symmetric and positive definite.
    r = b - A * x0;
    p = r;
    rsold = r' * r;
    x = x0;
    while sqrt(rsold) > tol
        Ap = A * p;
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
end