function [p_final, l_final, u_final] = lu_decomp_partial_pivot(a)

    % this function return the LU decomposition from partial pivoting.
    
    mat_dim = size(a);
    n = mat_dim(1);
    m = a;
    l_final = eye(mat_dim);
    p_final = eye(mat_dim);
    for i=1:n-1
        temp1 = max(m(i:n,i));
        temp2 = min(m(i:n,i));
        if abs(temp1) < abs(temp2)
            pos = find(m(i:n,i) == temp2) +(i-1);
        else
            pos = find(m(i:n,i) == temp1) +(i-1);
        end
        if pos~=i
            p = eye(mat_dim);
            p(i,i) = 0;
            p(i,pos) = 1;
            p(pos,pos) = 0;
            p(pos,i) = 1;
            p_final = p * p_final;
            for k = 1:i-1
                temp_1 = l_final(pos,k);
                l_final(pos,k) = l_final(i,k);
                l_final(i,k) = temp_1;
            end
            m = p * m;
        end
        % till here we computed permutation matrice and updated biggest
        % element as pivot.next we do gaussian elimination by explicitly
        % saving coeffcients of eliminatoin in L matrice at last we get U
        % matrice. and we just irirate this process to get whole matrice as
        % LU
        
        for j = n:-1:i+1
            e = eye(mat_dim);
            coff = m(j,i)/m(i,i);
            e(j,i) = - coff;
            m = e * m;
            l_final(j,i) = coff;
        end
    end
    u_final = m;
end
