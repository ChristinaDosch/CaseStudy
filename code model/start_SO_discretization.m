function [x_opt, obj_opt, runningTime] = start_SO_discretization
%% DISCRETIZATION APPROACH
% This script solves our well known optimization problem using the
% discretization approach and including
%       * ramping constraints and x_min, x_max
%       * no battery constraints
%       * general penalty function (also quadratic is possible)

%% Initialize parameters
[T, P, cost, penalty, epsilon, C] = init_parameters;

%% Constraints
[x_min,x_max,delta,A,b] = init_constraints(T,P,C);

%% Example scenarios

% for this example T = 1440 is required in init_parameters!!!
load('PVdata2'); 
K = 31;
E = reshape(PVdata2(:,1),31,1440);
%K = 372; % number of realizations
%E = reshape(PVdata2,372,1440); % array of K realizations (one per row) with data per minute of one day, respectively

F = cell(K,1); %F(k) contains F(x,e^k)

% determine revenue function F(x,e^k) for every e^1,...,e^K:
for k = 1:K
e = E(k,:); % contains the k^th realization
F(k) = {@(x) obj_SO_discr(x,e,cost,penalty,epsilon,P)}; 
end

objfct = @(x) 1/K * sum(cellfun(@(f)f(x),F)); % weighted (all weights=1/K) sum of F(x,e^k)

%% Performing optimization
x0 = zeros(1,T);
tic
[x_opt, obj_opt] = fmincon(objfct,x0,A,b,[],[],x_min*ones(1,T),x_max*ones(1,T)); 
runningTime = toc
end