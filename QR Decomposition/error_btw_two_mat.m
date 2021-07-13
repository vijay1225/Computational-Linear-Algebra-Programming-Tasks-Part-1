function e = error_btw_two_mat(a,b)

    % this function just returns the norm(A - B) vectors which is error
    
    mat_dim = size(a);
    sum = 0;
    for i=1:mat_dim(1)
        for j=1:mat_dim(2)
            sum = sum + (a(i,j)-b(i,j))^2;
        end
    end
    e = sqrt(sum);
end