clc;
clear;

% this code generate toeplitz matrice for given n and estimate AX=b with
% jacobi and gauss seidel methods
% i have used fuctions which return the sol of x and cound along with
% convergence rates

disp('since toepliz matrice is diagonally dominant it will surely be converge');
for n = [10 20]
    disp('for n =');
    disp(n);
    col = zeros(n,1);
    col(1,1) = 2;
    col(2,1) = -1;
    a = toeplitz(col);
    b = zeros(n,1);
    b(1,1) = 11;
    epsilon = 0.0001;

    [x_jacobi, count_jacobi,plot_jacobi, x_gauss, count_gauss, plot_gauss] = jacobi_gauss_1(a,b,epsilon);

    x = a\b;
    disp('Error estimate with x');
    disp('error to toepliz matrice to given n value with jacobi method');
    disp(norm(x - x_gauss));
    disp('error to toepliz matrice to given n value with gauss method');
    disp(norm(x - x_gauss));
end

