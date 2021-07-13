function [x_jacobi, count_jacobi, plot_jacobi_values, x_gauss, count_gauss, plot_gauss_values] = jacobi_gauss_1(a,b,epsilon)
    
    % this function returns the both jacobi and gauss seidel iterative
    % method and returns the iterations cound and convergence rate 
    % ||x^k - x^(k-1)|| <= (sigma)^k * ||x^1 - x^0|| is the stopping
    % critiran
    
    dim = size(a);
    n = dim(1);
    l = tril(a,-1);
    u = triu(a,1);
    d = a-l-u;
    plot_jacobi_values = [];
    plot_gauss_values = [];
    %% Jacobi Method
    
    % since diagonal matrice inverse is same as taking individual element
    % reciprocal i have used matrice multiplication update rule only
    
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
    
    
    
    q_gauss = (l+d)\b;
    p_gauss = - (l+d)\u;

    x_gauss = zeros(n,1);
    x_gauss_iter1 = (p_gauss * x_gauss) + q_gauss;
    plot_gauss_values(1) = sqrt((x_gauss_iter1 - x_gauss)' * (x_gauss_iter1 - x_gauss));
    tolarance = epsilon * plot_gauss_values(1);
    x_gauss = x_gauss_iter1;
    count = 1;
    while true
        x_gauss_old = x_gauss;
% 
%             for i=1:n
%                 sum = 0;
%                 for j=1:i-1
%                     sum = sum + (a(i,j) + x_gauss(j,1));
%                 end
%                 for j = i+1:n
%                     sum = sum + (a(i,j) + x_gauss_old(j,i));
%                 end
%                 x_gauss(i,1) = (-(sum)/a(i,i))+(b(i,1)/a(i,i));
%             end
%             count = count + 1;
%             if count == 1:
%                 plot_gauss_values(1) = sqrt((x_gauss - x_gauss_old)' * (x_gauss - x_gauss_old));
%                 tolarance = epsilon * plot_gauss_values(1);
%             end
        
        % trial
        
        x_gauss = (p_gauss * x_gauss) + q_gauss;
        count = count + 1;
        check = sqrt((x_gauss - x_gauss_old)' * (x_gauss - x_gauss_old));
        plot_gauss_values(count) = check;
        if  check < tolarance
            break
        end
    end
    count_gauss = count;
end