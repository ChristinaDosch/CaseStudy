function [x_min, x_max, delta, A, b, A_, b_, A_smart, b_smart] = init_constraints(T, P, C, SOC_0)
% returns bounds x_min,x_max,delta for constraints as well as matrix A and
% vector b, describing all constraints in the form Ax<=b.
% For detailled description of A and b see "Ansatz 2: Robust Optimization"
% 
% input T indicates the number of time steps
% input P indicates nominal power (needed since delta and x_max depend on
% P)
% input C: battery capacity

x_min = 0;  % minimum amount of power that can be scheduled
x_max = 0.7*P;  % maximum amount of power that can be scheduled
delta = 0.03*P;    % maximum deviation between two successive power values
B = [-eye(T-1) zeros(T-1,1)] + [zeros(T-1,1) eye(T-1)];
A = [B; -B];
b = ones(2*(T-1),1)*delta;

% For RO with battery
A_ = tril(ones(T)); A_ = [A_; -A_];
b_ = [C*ones(T,1); zeros(T,1)];

% For SO with battery (SMART)
B_tilde = [B, zeros(T-1,3*T)];
C_smart = [zeros(T,T), eye(T) + diag(-ones(T-1,1),-1) , -0.95/C * eye(T), 1/C * eye(T)];
D = [zeros(T,2*T), -eye(T), 0.95 * eye(T)];
c = [SOC_0; zeros(T-1,1)];
% d = ; %kann man hoffentlich weglassen
A_smart = [B_tilde; -B_tilde; C_smart; -C_smart; D; -D];
b_smart = [ones(T-1,1)*delta; -ones(T-1,1)*delta; c; -c; d; -d];

%TO DO:
% - upper bounds und lower bounds (die ja nicht in Matrix mit drin stehen,
%   weil Matlab mit solchen capacity constraints sinnvoller umgehen kann)
%   müssen dann auch oben als function output noch rein
% - d weglassen und damit auch D
end

