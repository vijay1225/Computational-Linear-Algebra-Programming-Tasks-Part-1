clc;
clear;

% this code prints the given system solutions with given tolerances.

% input to modify
epsilon = [0.1, 0.01, 0.001, 0.0001, 0.00001];
n = [5, 50, 100];

% code starts here
epsilon_size = size(epsilon);
n_size = size(n);
total_count_array = zeros(n_size(2), epsilon_size(2));
% iterated throuhg each and and each epsilon values
for i=1:n_size(2)
    count_array = zeros(1,epsilon_size(2)); 
    e=ones(n(1,i),1);
    A=spdiags([-e 2*e -e], -1:1, n(1,i), n(1,i));
    A=full(A);
    b=rand(n(1,i),1);
    for j=1:epsilon_size(2)
        [x_jacobi, count, time] = jacobi_method(A,b,epsilon(1,j));
        count_array(1,j) = count;
    end
    total_count_array(i,:) = count_array;
    plot(1:5, log(count_array));
    title('Tolarance vs count');
    xlabel('Tolarance (log scale)');
    ylabel('count(log scale)');
    hold on
end
disp(total_count_array);
hold off

disp('All plots are taken in logorithmic scale in both axis to get good understanding');
disp(' as n increasing itirations needed to get good accuracy increases we can observe it from plot');
