function [x_min, x_max, delta, A, b, C, A_, b_] = init_constraints(T)
% returns bounds x_min,x_max,delta for constraints as well as matrix A and
% vector b, describing all constraints in the form Ax<=b.
% For detailled description of A and b see "Ansatz 2: Robust Optimization"
% input T indicates the number of time steps

x_min = 0;  % minimum amount of power that can be scheduled
x_max = 2.25;  % maximum amount of power that can be scheduled
delta = 0.15;    % maximum deviation between two successive power values
B = [-eye(T-1) zeros(T-1,1)] + [zeros(T-1,1) eye(T-1)];
A = [B; -B];
b = ones(2*(T-1),1)*delta;

% F?r RO with battery
C = 3; % capacity
A_ = tril(ones(T)); A_ = [A_; -A_];
b_ = [C*ones(T,1); zeros(T,1)];
end

