function [q,r] = hr(input_mat)
    
    % i have created a reflections functoin that retuns the sub matice
    % after doing operation of householder transformation.
    
    % took input matrice and reduce it into uppertriangular by doing
    % blockwise householder transformantion.

    mat_dim = size(input_mat);
    n = mat_dim(2);
    r = input_mat;
    q = eye(mat_dim);
    for i=1:n
       t = eye(mat_dim);
       t(i:n,i:n) = reflection(r(i:n,i:n));
       q = q * t';
       r = q' * input_mat;
    end
end
