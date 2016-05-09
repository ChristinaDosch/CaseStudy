function [x_min, x_max, delta, A, b] = constraints(N)
% returns bounds x_min,x_max,delta for constraints as well as matrix A and
% vector b, describing all constraints in the form Ax<=b.
% For detailled description of A and b see "Ansatz 2: Robust Optimization"
% input N indicates the number of time steps

x_min = 0;  % minimum amount of power that can be scheduled
x_max = 0.55; % maximum amount of power that can be scheduled
delta = 0.04; % maximum deviation between two successive power values
B = [-eye(N-1) zeros(N-1,1)] + [zeros(N-1,1) eye(N-1)];
A = [B; -B];
b = ones(2*(N-1),1)*delta;

end

