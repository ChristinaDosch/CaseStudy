function [x_min, x_max, delta, SOC_min, SOC_max, A, b, A_, b_, A__, b__] = init_constraints_inout(T,...
    P, C, SOC_0)
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
%          A_,b_: contain constraints for b_out with battery
%          A__,b__: contain capacity constrains

% capacity bounds:
x_min = 0;      % minimum amount of power that can be scheduled
x_max = 0.7*P;  % maximum amount of power that can be scheduled
SOC_min = 0.1;  % maximum of state of charge
SOC_max = 0.95; % minimum of state of charge

% Ramping constraints for x (not capacity constraints)
% delta = 0.045 * P; % maximum deviation between two successive power values in case of T=96
delta = 0.02*P;
% delta = 0.01*P;

if T == 25, delta = 4*delta; end % ... in case of T=25

B = [-eye(T-1) zeros(T-1,1)] + [zeros(T-1,1) eye(T-1)];
A = [B zeros(T-1,T) zeros(T-1,T);
    -B zeros(T-1,T) zeros(T-1,T)];
b = ones(2*(T-1),1)*delta;

% Constraints for b_out: bi_out <= b0 + sum[k=1..i-1](0.95*bk_in - bk_out)
A_ = [zeros(T) -0.95*(tril(ones(T))-eye(T)) tril(ones(T))];
b_ = ones(T,1)*SOC_0*C;


% Capacity constraints for SOC*C for RO
A__ = [zeros(T)  0.95*tril(ones(T)) -tril(ones(T));
       zeros(T) -0.95*tril(ones(T))  tril(ones(T))];
b__ = [ones(T,1)*(SOC_max - SOC_0)*C;
       ones(T,1)*(SOC_0 - SOC_min)*C];

end