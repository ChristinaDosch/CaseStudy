function [x_min, x_max, SOC_min, SOC_max, delta, A, b, A_, b_, A_smart, b_smart] = init_constraints(T, P, C, SOC_0)
% returns bounds x_min, x_max, delta for constraints as well as matrix A and
% vector b, describing all constraints in the form Ax<=b.
% For detailled description of A and b see "Ansatz 2: Robust Optimization"
% 
% input:   T - scalar:           number of time steps
%          P - scalar:           nominal power (needed since delta and x_max depend on P)
%          C - scalar:           battery capacity

% capacity bounds:
x_min = 0;      % minimum amount of power that can be scheduled
x_max = 0.7*P;  % maximum amount of power that can be scheduled
SOC_min = 0.1;  % maximum of state of charge
SOC_max = 0.95; % minimum of state of charge

% "real" constraints (not capacity constraints)
delta = 0.03*P; % maximum deviation between two successive power values
B = [-eye(T-1) zeros(T-1,1)] + [zeros(T-1,1) eye(T-1)];
A = [B; -B];
b = ones(2*(T-1),1)*delta;

% For RO with battery
A_ = tril(ones(T)); A_ = [A_; -A_];
b_ = [C*ones(T,1); zeros(T,1)];

% For SO with battery (SMART)
B_tilde = [B, zeros(T-1,3*T)];
C_smart = [zeros(T,T), eye(T) + diag(-ones(T-1,1),-1) , -0.95/C * eye(T), 1/C * eye(T)];
c = [SOC_0; zeros(T-1,1)];

A_smart = [B_tilde; -B_tilde; C_smart; -C_smart];
b_smart = [ones(T-1,1)*delta; ones(T-1,1)*delta; c; -c];

%TO DO:
% - upper bounds und lower bounds (die ja nicht in Matrix mit drin stehen,
%   weil Matlab mit solchen capacity constraints sinnvoller umgehen kann)
%   müssen dann auch oben als function output noch rein -erledigt
% - d weglassen und damit auch D -erledigt
end

