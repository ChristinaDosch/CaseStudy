function [x_min, x_max, delta, A, b, A_, b_, A_smart, b_smart, SOC_min, SOC_max] = init_constraints(T, P, C, SOC_0,K)
%% INIT_CONSTRAINTS(T,P,C,SOC_0) returns all constraints in the form Ax<=b
% 
% input:   T - scalar:                  number of time steps
%          P - scalar:                  nominal power (needed since delta and x_max depend on P)
%          C - scalar:                  battery capacity
%          SOC_0 - scalar in [0,1]:     state of charge at day break
%
%          all initialized in init_parameters
%
% output:  upper bounds x_max and SOC_max and lower bounds x_min and SOC_min for variables x and SOC, respectively
%          delta: parameter of ramping constraint
%          A, b: contain constraints for model without battery (described by Ax<=b)
%          A_,b_: contain constraints for RO with battery
%          A_smart, b_smart: contain constraints for SO with smart battery
%
%          for detailed derivation of these matrices and vectors see
%          chapter 2.2 in documentation

% capacity bounds:
x_min = 0;      % minimum amount of power that can be scheduled
x_max = 0.7*P;  % maximum amount of power that can be scheduled
SOC_min = 0.1;  % maximum of state of charge
SOC_max = 0.95; % minimum of state of charge

% Ramping constraints for x (not capacity constraints)
delta = 0.03 * P; % maximum deviation between two successive power values
B = [-eye(T-1) zeros(T-1,1)] + [zeros(T-1,1) eye(T-1)];
A = [B; -B];
b = ones(2*(T-1),1)*delta;

% Capacity constraints for SOC*C
A_ = tril(ones(T)); A_ = [A_; -A_];
b_ = [C*ones(T,1); zeros(T,1)];

% For SO with battery (SMART), requires declaration of K in the function call
if nargin == 5
    % matrix for ramping constraints:    
    B_tilde = [B, zeros(T-1,3*K*T)];
    % Nebendiagonale for submatrix of C:
    nd = -ones(T*K-1,1); 
        for i=1:K-1, nd(i*T) = 0; end 
    % matrix for SOC-constraints:    
    C_smart = [zeros(K*T,T), eye(K*T) + diag(nd,-1) , -0.95/C * eye(K*T), 1/C * eye(K*T)];
    % right hand vector for SOC-constraints:
    c = zeros(K*T,1);
        for i=0:K-1, c(i*T+1) = SOC_0; end
    
    % entire matrix-vector system for smart battery with discretization    
    A_smart = [B_tilde; -B_tilde; C_smart; -C_smart];
    b_smart = [ones(T-1,1)*delta; ones(T-1,1)*delta; c; -c];
end
end