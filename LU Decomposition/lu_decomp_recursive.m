function [l_final,u_final] = lu_decomp_recursive(a)
    
    % this function return the LU decomposition using recursive method
    % which is based on separating L and U in the form of Rank 1
    % matriaces.
    
    % we separate L and U columns and row each, and find block of which we again
    % find L U colums and rows iteratively we get LU matrices.
    
    mat_dim = size(a);
    n = mat_dim(1);
    m = a;
    l_final = zeros(mat_dim);
    u_final = zeros(mat_dim);
    for i=1:n
        l = m(:,i);
        u = m(i,:);
        if abs(l(i,1)) < 0.000001
            disp('LU decomposition not exist');
            break;
        end
        l = l/l(i,1);
        u_final(i,:) = u;
        l_final(:,i) = l;
        m = m - (l * u);
    end
end
    
    
    
    