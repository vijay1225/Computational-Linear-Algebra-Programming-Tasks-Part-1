clc;
clear;

epsilon = 0.001;
n = [5, 50, 100];

n_size = size(n);
total_time_array = zeros(n_size(2),1);
for i = n
    e=ones(i,1);
    A=spdiags([-e 2*e -e], -1:1, i, i);
    A=full(A);
    b=rand(i,1);
    [x_jacobi, count_jacobi, time_jacobi] = jacobi_method(A,b,epsilon);
    tic;
    sol = cmplx(A,b);
    time_direct_method = toc;
    disp('for n =');
    disp(i);
    disp('time taken to process with jacobi and LU decompostion respectively');
    disp(time_jacobi);
    disp(time_direct_method);
end
disp('we can conclude that instead of iteratiive method direct method is faster.');
disp('but obtaining direct method matrices always is may not possible');

