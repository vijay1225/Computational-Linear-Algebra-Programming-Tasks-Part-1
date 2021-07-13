function [sol] = cmplx(a,b)

    % return the system solution using lu decompositon and back
    % substitution
    % lu_decom_recursice funtionis used which decompose A into LU
    % recursively
    
    b_size = size(b);
    if b_size(1) == 1
        b = b';
        b_size = size(b);        
    end
    [l,u] = lu_decomp_recursive(a);
    % forward substitution
    y = zeros(b_size);
    for i=1:b_size(1)
        sum = 0;
        for j = 1:i-1
            sum = sum + l(i,j) * y(j,1);
        end
        y(i,1) = (b(i,1) - sum)/l(i,i);
    end
    % backward substitution
    x = zeros(b_size);
    for i=b_size(1):-1:1
        sum = 0;
        for j = b_size(1):-1:i+1
            sum = sum + u(i,j) * x(j,1);
        end
        x(i,1) = (y(i,1) - sum)/u(i,i);
    end
    sol = x;
end
