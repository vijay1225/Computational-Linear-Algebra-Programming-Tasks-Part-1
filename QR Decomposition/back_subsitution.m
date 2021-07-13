function [sol] = back_subsitution(b,u)

    % this function returns the solution of Ux = b by back substitution
    % since U is upper triangular matrice.(only works if U is upper
    % triangular
    
    mat_dim = size(b);
    n = mat_dim(1);
    x = zeros(n,1);
    for i=n:-1:1
        if u(i,i) > 0.0000001
            sum = 0;
            for j =n:-1:i+1
                sum = sum + u(i,j) * x(j,1);
            end

            x(i,1) = (b(i,1) - sum)/u(i,i);
        end
    end
    sol = x;
end
