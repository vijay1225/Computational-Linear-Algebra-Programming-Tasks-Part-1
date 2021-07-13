function [x_jacobi, count_jacobi, plot_jacobi_values, x_gauss, count_gauss, plot_gauss_values] = jacobi_gauss_1(a,b,epsilon)
    dim = size(a);
    n = dim(1);
    l = tril(a,-1);
    u = triu(a,1);
    d = a-l-u;
    plot_jacobi_values = [];
    plot_gauss_values = [];
    %% Jacobi Method
    % As in jacobi method D is diagonal its inverse is computational
    % efficient so we can use direct x = p * x + b formula
    q_jacobi = d\b;
    p_jacobi = - d\(l+u);

    x_jacobi = zeros(n,1);
    x_jacobi_iter1 = (p_jacobi * x_jacobi) + q_jacobi;
    plot_jacobi_values(1) = sqrt((x_jacobi_iter1 - x_jacobi)' * (x_jacobi_iter1 - x_jacobi));
    tolarance = epsilon * plot_jacobi_values(1);
    x_jacobi = x_jacobi_iter1;
    count = 1;
    while true
        x_jacobi_old = x_jacobi;
        x_jacobi = (p_jacobi * x_jacobi) + q_jacobi;
        count = count + 1;
        check = sqrt((x_jacobi - x_jacobi_old)' * (x_jacobi - x_jacobi_old));
        plot_jacobi_values(count) = check;
        if  check < tolarance
            break
        end
    end
    count_jacobi = count;
    
    %% Guass Seidel Method
    % As gauss seidal method needs inv(L+U) i followed element wise
    % updation
    x_gauss = zeros(n,1);
    tolarance = epsilon * plot_gauss_values(1);
    count = 0;
    while true
        x_gauss_old = x_gauss;
        for i=1:n
            sum = 0;
            for j=1:i-1
                sum = sum + (a(i,j) + x_gauss(j,1));
            end
            for j = i+1:n
                sum = sum + (a(i,j) + x_gauss_old(j,i));
            end
            x_gauss(i,1) = (-(sum)/a(i,i))+(b(i,1)/a(i,i));
        end
        count = count + 1;
        if count == 1:
            plot_gauss_values(1) = sqrt((x_gauss - x_gauss_old)' * (x_gauss - x_gauss_old));
            tolarance = epsilon * plot_gauss_values(1);
        end
        check = sqrt((x_gauss - x_gauss_old)' * (x_gauss - x_gauss_old));
        plot_jacobi_values(count) = check;
        if  check < tolarance
            break
        end
    end
