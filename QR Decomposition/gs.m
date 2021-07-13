function [q,r] = gs(input_mat)

    %taking inputs and preallocating orthobasis to zeros
    mat_dim = size(input_mat);
    ortho_mat = zeros(mat_dim);
    % following is algorithm is implemented
    % 1)first for loop takes the input matrice colums one by one as v
    % 2)second for loop takes already found orthovectors and find projection
    % wih selected column and sum it up
    % 3)now new vector is normalized and added to othobasis.
    for j = 1:mat_dim(2)
        v = transpose(input_mat(:,j));
        e = zeros(mat_dim(1),1);
        for l = 1:j-1
            e = e + (v * ortho_mat(:,l) * ortho_mat(:,l));
        end
        e = transpose(v) - e;
        if transpose(e) * e > 2 * 10^-16
            e = e/sqrt(transpose(e)* e);
            ortho_mat(:,j) = e;
        end
    end
    r = transpose(ortho_mat) * input_mat;
    q = ortho_mat;
end