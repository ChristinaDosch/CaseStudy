function [x_opt, obj_opt, runningTime] = start_SO_discretization_battery
%% DISCRETIZATION APPROACH with BATTERY
% This script solves our well known optimization problem using the
% discretization approach and including
%       * ramping constraints and x_min, x_max
%       * battery using the brute force approach
%       * general penalty function (also quadratic is possible)

%% Initialize parameters
[T, P, cost, penalty, epsilon, C] = init_parameters;

%% Constraints
[x_min,x_max,delta,A,b] = init_constraints(T,P,C);

%% Example scenarios
load('PVdata2');
%K = 372; % number of realizations
%E = reshape(PVdata2,372,1440); % array of K realizations (one per row) with data per minute of one day, respectively
K = 31;
E = reshape(PVdata2(:,1),31,1440); % array of 31 realizations with minute values from January
F = cell(K,1);

% determine revenue function F(x,\tilde{x}^k) for every k=1,...,K:
for k = 1:K 
x_tilde = @(x) battery(E(k,:),x,C); % compute \tilde{x}^k
F(k) = { @(x) obj_SO_discr(x,x_tilde(x),cost,penalty,epsilon,P)};
end

objfct = @(x) 1/K * sum(cellfun(@(f)f(x),F)); % weighted (all weights=1/K) sum of F(x,e^k)

%% Performing optimization
x0 = zeros(1,T);
tic
[x_opt, obj_opt] = fmincon(objfct,x0,A,b,[],[],x_min*ones(1,T),x_max*ones(1,T));
runningTime = toc
end