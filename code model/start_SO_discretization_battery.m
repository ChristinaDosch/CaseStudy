%% DISCRETIZATION APPROACH
% This script solves our well known optimization problem using the
% discretization approach and including
%       * ramping constraints and x_min, x_max
%       * battery using the brute force approach
%       * general penalty function (also quadratic is possible)

%% Initialize parameters
[T, P, cost, penalty, epsilon, t, mu, sigma, lambda] = init_parameters;
C = 2.6; % battery capacity (to be included in init_parameters)

%% Constraints
x_min = 0;  
x_max = 0.7*P; 
delta = 0.03*P; % maximum deviation allowed
B = [-eye(T-1) zeros(T-1,1)] + [zeros(T-1,1) eye(T-1)];
A = [B; -B];
b = ones(2*(T-1),1)*delta;

%% Example scenarios
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
tic
[x_opt, obj_opt] = fmincon(objfct,x0,A,b,[],[],x_min*ones(1,T),x_max*ones(1,T));
toc
