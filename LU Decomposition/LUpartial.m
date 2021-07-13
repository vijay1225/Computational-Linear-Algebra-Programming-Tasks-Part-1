function [a_inv, determinant] = LUpartial(p,l,u)

    % this function returns the inverse(A) to given LUP matrices ny
    % considering each colum of inverse is a estimate LUx=P
    % to get whole inverse we iterate through all colums of inverse and
    % estimated result

    determinant = det(l) * det(u);
    mat_dim = size(p);
    if determinant ~= 0
        n = mat_dim(1);
        a_inv = zeros(mat_dim);
        for t=1:n
            b = p(:,t);
            y = zeros(n,1);
            % forward substitution
            for i=1:n
                sum = 0;
                for j = 1:i-1
                    sum = sum + l(i,j) * y(j,1);
                end
                y(i,1) = (b(i,1) - sum)/l(i,i);
            end
            % backward substitution
            x = zeros(n,1);
            for i=n:-1:1
                sum = 0;
                for j =n:-1:i+1
                    sum = sum + u(i,j) * x(j,1);
                end
                x(i,1) = (y(i,1) - sum)/u(i,i);
            end
            a_inv(:,t) = x;
        end
    else
        disp('matrice not invertible');
        a_inv = zeros(mat_dim);
    end
end