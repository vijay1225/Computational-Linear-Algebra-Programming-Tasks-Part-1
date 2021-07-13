function [x_jacobi, count_jacobi,time] = jacobi_method(a,b,epsilon)

    % this function perform the jacobi method to given system and returns
    % time taken to process and no of iterations.
    tic;
    dim = size(a);
    n = dim(1);
    l = tril(a,-1);
    u = triu(a,1);
    d = a-l-u;
    %% Jacobi Method
    q_jacobi = d\b;
    p_jacobi = - d\(l+u);

    x_jacobi = zeros(n,1);
    x_jacobi_iter1 = (p_jacobi * x_jacobi) + q_jacobi;
    check = sqrt((x_jacobi_iter1 - x_jacobi)' * (x_jacobi_iter1 - x_jacobi));
    tolarance = epsilon * check;
    x_jacobi = x_jacobi_iter1;
    count = 1;
    while true
        x_jacobi_old = x_jacobi;
        x_jacobi = (p_jacobi * x_jacobi) + q_jacobi;
        count = count + 1;
        check = sqrt((x_jacobi - x_jacobi_old)' * (x_jacobi - x_jacobi_old));
        if  check < tolarance
            break
        end
    end
    count_jacobi = count;
    time = toc;
end